import time
import os
import multiprocessing
import concurrent.futures
import argparse

from tstools.ts_optimizer import TSOptimizer
from tstools.utils import remove_files_in_directory, copy_final_outputs, \
    setup_dir, get_reaction_list, print_statistics


def get_args():
    """
    Parse command-line arguments.

    Returns:
    - argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--reactive-complex-factors-intra', nargs='+', type=float,
                        default=[0, 1.2, 1.3, 1.8])
    parser.add_argument('--reactive-complex-factors-inter', nargs='+', type=float, 
                        default=[2.5, 1.8, 2.8, 1.3])
    parser.add_argument('--freq-cut-off', action='store', type=int, default=150)
    parser.add_argument('--solvent', action='store', type=str, default=None)
    parser.add_argument('--xtb-external-path', action='store', type=str, 
                        default="xtb_external_script/xtb_external.py")
    parser.add_argument('--input-file', action='store', type=str, default='data/reactions_am.txt')
    parser.add_argument('--target-dir', action='store', type=str, default='work_dir')

    return parser.parse_args()


def optimize_individual_ts(ts_optimizer):
    """
    Optimize an individual transition state.

    Parameters:
    - ts_optimizer: Instance of TSOptimizer.

    Returns:
    - int or None: Reaction ID if a transition state is found, None otherwise.
    """
    # First select the set of reactive_complex factor values to try
    start_time_process = time.time()

    try:    
        if ts_optimizer.reaction_is_intramolecular():
            reactive_complex_factor_values = ts_optimizer.reactive_complex_factor_values_intra
        else:
            reactive_complex_factor_values = ts_optimizer.reactive_complex_factor_values_inter
    except Exception as e:
        print(e)
        return None

    # Then search for TS by iterating through reactive complex factor values
    for reactive_complex_factor in reactive_complex_factor_values:
        for _ in range(3):
            try:
                ts_optimizer.set_ts_guess_list(reactive_complex_factor)
                ts_found = ts_optimizer.determine_ts() 
                remove_files_in_directory(os.getcwd())
                if ts_found:
                    end_time_process = time.time()
                    print(f'Final TS guess found for {ts_optimizer.rxn_id} for reactive complex factor {reactive_complex_factor} in {end_time_process - start_time_process} sec...')
                    return ts_optimizer.rxn_id
            except Exception as e:
                print(e)
                continue

    end_time_process = time.time()
    print(f'No TS guess found for {ts_optimizer.rxn_id}; process lasted for {end_time_process - start_time_process} sec...')

    return None


def obtain_transition_states(target_dir, reaction_list, xtb_external_path, solvent,
                             reactive_complex_factor_list_intermolecular,
                             reactive_complex_factor_list_intramolecular, freq_cut_off):
    """
    Obtain transition states for a list of reactions.

    Parameters:
    - target_dir (str): Target directory.
    - reaction_list (list): List of reactions.
    - xtb_external_path (str): Path to the XTB external script.
    - solvent (str): Solvent information.
    - reactive_complex_factor_list_intermolecular (list): List of reactive complex factors for intermolecular reactions.
    - reactive_complex_factor_list_intramolecular (list): List of reactive complex factors for intramolecular reactions.
    - freq_cut_off (int): Frequency cutoff.

    Returns:
    - list: List of successful reactions.
    """
    home_dir = os.getcwd()
    os.chdir(target_dir)
    ts_optimizer_list = []

    for rxn_idx, rxn_smiles in reaction_list:
        ts_optimizer_list.append(TSOptimizer(rxn_idx, rxn_smiles, xtb_external_path,
                                             solvent, None, reactive_complex_factor_list_intermolecular,
                                             reactive_complex_factor_list_intramolecular, freq_cut_off))

    print(f'{len(ts_optimizer_list)} reactions to process...')

    num_processes = multiprocessing.cpu_count()

    with concurrent.futures.ProcessPoolExecutor(max_workers=int(num_processes/2)) as executor:
        # Map the function to each object in parallel
        results = list(executor.map(optimize_individual_ts, ts_optimizer_list))

    successful_reactions = [r for r in results if r is not None]

    os.chdir(home_dir)
    copy_final_outputs(target_dir, f'final_{target_dir}')

    return successful_reactions
    

if __name__ == "__main__":
    # preliminaries
    args = get_args()
    setup_dir(args.target_dir)
    reaction_list = get_reaction_list(args.input_file)
    start_time = time.time()
    xtb_external_path = f'{os.path.join(os.getcwd(), args.xtb_external_path)}'

    # run all reactions in parallel
    successful_reactions = obtain_transition_states(args.target_dir, reaction_list, 
        xtb_external_path, solvent=args.solvent, 
        reactive_complex_factor_list_intramolecular=args.reactive_complex_factors_intra, 
        reactive_complex_factor_list_intermolecular=args.reactive_complex_factors_inter, 
        freq_cut_off=args.freq_cut_off)
    
    # print final statistics about the run
    print_statistics(successful_reactions, start_time)
