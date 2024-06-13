#!/bin/bash
#SBATCH -J Strychnine
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH --export=NONE
#SBATCH --ntasks-per-node=24

##################################################
################    MODULES     ##################
##################################################

source $HOME/.bash_profile

g16root="/opt/apps/Gaussian/G16C"
export g16root
. $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR="/tmp"

source activate DL_CPU

##################################################
################    VARIABLES    #################
##################################################

childcpu=12 # Number of cpus per parallel job during batching
totalcpu=$SLURM_NTASKS_PER_NODE
jbatch=$(expr $totalcpu / $childcpu)
excdir=$(pwd)

##################################################
################    FUNCTIONS    #################
##################################################

### Run QM calculation. In theory any program could be called here.
task () {

$g16root/g16/g16 < $1.com >> $2.log

}

### Wrapper for batching  and running jobs
gexc () {
#Variables; Use local here or otherwise it becomes global and overwrites original
local jobname=$1 # The current job
local jbatch=$2 # Number or parallel jobs.

if [ -d $jobname ]; then
        cd $jobname # We don't exit this directory within the function. Might cause issues.
else
        echo "No $jobname folder found, exiting the script"
        exit
fi
start=$(date +%s)
for i in $(ls -1vd *.com); do
        (
        file=$(echo $i | sed 's/.com//')
        echo "$file"
        task $file $excdir/$jobname/$file
        ) &
        # allow to execute up to $N jobs in parallel
        if [[ $(jobs -r -p | wc -l) -ge $jbatch ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
                wait
        fi
done
wait
end=$(date +%s)
sec=$(expr $end - $start)
echo "Time $jobname: $sec seconds"

}


##################################################
###############    MAIN SCRIPT     ###############
##################################################


# Run CSEARCH conformational search
python -m aqme --csearch --program rdkit --smi 'C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75' --name Strychnine

# Run QPREP to create Gaussian input files
python -m aqme --qprep --program gaussian --files "CSEARCH/*sdf" --qm_input 'B3LYP/6-31+G(d,p) opt freq' --mem '24GB' --nprocs "$childcpu"

# Run Gaussian with the generated input files
echo "Starting opt-freq calculations"
gexc QCALC $jbatch
wait
cd $excdir

# Run first QCORR analysis
python -m aqme --qcorr --files "$excdir/QCALC/*log" --mem '24GB' --nprocs "$childcpu" --freq_conv 'opt=(calcfc,maxstep=5)' --isom_type 'com' --isom_inputs "$excdir/QCALC"

# Run second round of Gaussian/QCORR analysis
if [ -d "$excdir/QCALC/failed/run_1/fixed_QM_inputs" ]; then
    echo "Starting one rerun of failed calculations"
            gexc $excdir/$jobname/$fixed $jbatch
            wait
            cd $excdir
            python -m aqme --qcorr --files "$excdir/QCALC/failed/run_1/fixed_QM_inputs/*.log" --fullcheck --mem '24GB' --nprocs "$childcpu" --isom_type 'com' --isom_inputs "$excdir/QCALC/failed/run_1/fixed_QM_inputs/*.log"
    else
            echo "No fixable errors found"
    fi

echo "Done with opt-freq calculations"
cd $excdir

#Create NMR input files to obtain chemical shift data
python -m aqme --qprep --program gaussian --files "$excdir/QCALC/success/*.log" --destination "$excdir/NMR-SP" --qm_input 'B3LYP/6-311+G(2d,p) scrf=(solvent=chloroform,smd) nmr=giao' --mem '24GB' --nprocs "$childcpu" --suffix 'SP'

# Run Gaussian with the generated NMR input files
echo "Starting SP calculations"
gexc NMR-SP $jbatch
wait
cd $excdir

#Run qcorr to create JSON files
python -m aqme --qcorr --program gaussian --files "$excdir/NMR-SP/*.log"

#Run Boltzmann averging of NMR with QDESCP
python -m aqme --qdescp --program nmr --boltz --files "$excdir/NMR-SP/success/SP_calcs/json_files/*_SP.json" --nmr_slope "[-1.0537, -1.0784]" --nmr_intercept "[181.7815,31.8723]"

# Run GoodVibes to obtain thermochemistry data, Boltzmann population and XYZ coordinates
python -m goodvibes --boltz --xyz -c 1 "$excdir/QCALC/success/*.log"
