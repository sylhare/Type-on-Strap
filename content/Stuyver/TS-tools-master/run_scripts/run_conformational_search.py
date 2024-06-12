import os
import shutil
import argparse
import subprocess
import numpy as np
from tstools.utils import write_final_geometry_to_xyz, extract_geom_from_xyz, extract_geom_from_crest_ensemble, create_input_file_opt_g16, extract_g16_energy


def get_args():
    """
    Parse command line arguments.

    Returns:
    -argparse.Namespace: Parsed command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, help='File to extract geometry')
    parser.add_argument('--opt-DFT', default=False, action="store_true", help='Subsequently optimization of the 10 conformers with the lowest energy, obtained from CREST to determine the most stable conformer')
    parser.add_argument('--atomlist', default=None, nargs="+", help='Constrained atoms')
    parser.add_argument('--charge', default=0, type=int, help='Charge of the molecular system')
    parser.add_argument('--uhf', default=0, type=int, help='Number of unpaired electrons of the molecular system')
    parser.add_argument('--mem', action='store', type=str, default='16GB', help='Specifies the memory requested in the Gaussian16 .com files')
    parser.add_argument('--proc', action='store', type=int, default=8, help='Number of CPU that will be used')    
    parser.add_argument('--functional', action='store', type=str, default='UB3LYP', help='Functional')
    parser.add_argument('--basis-set', action='store', type=str, default='6-31G', help='Basis set')
    parser.add_argument('--solvent', action='store', type=str, default=None, help='Solvent')    

    return parser.parse_args()


def run_crest(directory, name, geom, charge, uhf, atomlist, proc, solvent):
    """
    Run a conformational search using CREST.

    Parameters:
    -directory (str): Path to the directory for the calculation.
    -name (str): Name of the output file.
    -geom (list): List with the coordinates of the system.
    -charge (int): Charge of the molecular system.
    -uhf (int): Number of unpaired electrons of the molecular system.
    -atomlist (list): Constrained atoms
    -proc (int): Number of CPU that will be used 
    -solvent (str): Solvent

    Returns:
    None
    """    

    pwd = os.getcwd()
    os.chdir(directory)

    if solvent:
        kwd_solvent = f"--alpb {solvent}"
    else:
        kwd_solvent = " "

    if atomlist:
        print_constraints_inp(geom, atomlist)
        command_line = f"crest {name}.xyz --gfn2 -T {proc} --noreftopo --cinp constraints.inp {kwd_solvent} --chrg {charge} --uhf {uhf} > {name}.out"
    else:
        command_line = f"crest {name}.xyz --gfn2 -T {proc} --noreftopo {kwd_solvent}  --chrg {charge} --uhf {uhf} > {name}.out"

    with open('crest.out', 'w') as out:
        subprocess.run(command_line, shell=True, stdout=out, stderr=subprocess.DEVNULL,)
    os.chdir(pwd)


def print_constraints_inp(geom, atomlist):
    """
    Print the constraint file for CREST calculation.
    
    Parameters:
    -geom (list): List with the coordinates of the system.
    -atomlist (list): Constrained atoms

    Returns:
    None
    """
    constraints = []

    for i in range(0, len(atomlist), 2):
        idx_atm_1 = int(atomlist[i])
        idx_atm_2 = int(atomlist[i + 1])
        _, x1, y1, z1 = geom[idx_atm_1 - 1].split()
        _, x2, y2, z2 = geom[idx_atm_2 - 1].split()
        coord_atm_1 = np.array([x1, y1, z1], dtype=float)
        coord_atm_2 = np.array([x2, y2, z2], dtype=float)
        distance = np.linalg.norm(coord_atm_1 - coord_atm_2).round(3)
        constraints.append((idx_atm_1, idx_atm_2, distance))

    with open('constraints.inp', 'w') as file:
        file.write("$constrain\n")
        file.write("  force constant=0.5\n")
        for constraint in constraints:
            file.write(f"  distance: {constraint[0]}, {constraint[1]}, {constraint[2]}\n")
        file.write("$end")


def run_g16_opt(directory, name, charge, uhf, atomlist, mem, proc, functional, basis_set, solvent):
    """
     Run a conformational search using CREST.
 
    Parameters:
    -directory (str): Path to the directory for the calculation.
    -name (str): Name of the output file.
    -charge (int): Charge of the molecular system.
    -uhf (int): Number of unpaired electrons of the molecular system.
    -atomlist (list): Constrained atoms
    -mem (str): Specifies the memory requested in the Gaussian16 .com files
    -proc (int): Number of CPU that will be used
    -functional (str): functional
    -basis_set (str): basis_set
    -solvent (str): solvent
 
    Returns:
    None
    """
    
    os.chdir(new_dir)
    if atomlist:
        modredundant = formatting_constraints(atomlist)
    else:
        modredundant = None
    geoms = extract_geom_from_crest_ensemble('crest_conformers.xyz', 20)
    os.mkdir(f'{name}_opt_DFT')
    os.chdir(f'{name}_opt_DFT')
    multiplicity = int(2 * (uhf * 1/2) + 1)
    energies = []
 
    for idx, geom in enumerate(geoms):
        
        if not geom:
            break    
        tmp_name = f"{name}_conf_{idx}"
        command_line = f"g16 < {tmp_name}.com",
        create_input_file_opt_g16(name=tmp_name, geom=geom, charge=charge, multiplicity=multiplicity, mem=mem, proc=proc, modredundant=modredundant, functional=functional, basis_set=basis_set, solvent=solvent)
        with open(f'{tmp_name}.out' , 'w') as out:
            subprocess.run(command_line, shell=True, stdout=out, stderr=subprocess.DEVNULL,)

        energy = extract_g16_energy(f'{tmp_name}.out')
        
        if energy:
            energies.append((idx, energy))
    
    if energies:
        sorted_energies = sorted(energies, key=lambda x: x[1])
        os.mkdir('lowest')
        shutil.copy(f'{name}_conf_{sorted_energies[0][0]}.out', 'lowest')
    

def formatting_constraints(atomlist):

    constraints = []
    for i in range(0, len(atomlist), 2):
        idx_atm_1 = int(atomlist[i])
        idx_atm_2 = int(atomlist[i + 1])
        constraints.append(f"B  {idx_atm_1}  {idx_atm_2}  F")

    return constraints


if __name__ == '__main__':
    args = get_args()
    filename = args.file
    geom = extract_geom_from_xyz(filename)
    directory, name = os.path.split(filename)
    name = name.split('.')[0]
    new_dir = os.path.join(directory, f'{name}_crest')
    os.mkdir(new_dir)
    shutil.copy(filename, new_dir)
    run_crest(new_dir, name, geom, args.charge, args.uhf, args.atomlist, args.proc, args.solvent)    
    
    if args.opt_DFT:
        run_g16_opt(new_dir, name, args.charge, args.uhf,args.atomlist, args.mem, args.proc, args.functional, args.basis_set, args.solvent)
        


