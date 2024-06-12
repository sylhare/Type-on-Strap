import numpy as np
import os
import shutil
from rdkit import Chem

from tstools.path_generator import PathGenerator
from tstools.confirm_ts_guess import validate_ts_guess
from tstools.utils import xyz_to_gaussian_input, run_g16_ts_optimization, run_irc, remove_files_in_directory
from tstools.irc_search import generate_gaussian_irc_input, extract_transition_state_geometry, extract_irc_geometries, compare_molecules_irc

ps = Chem.SmilesParserParams()
ps.removeHs = False

class TSOptimizer:
    def __init__(self, rxn_id, reaction_smiles, xtb_external_path, xtb_solvent=None, dft_solvent=None,
                 reactive_complex_factor_values_inter=[2.5], reactive_complex_factor_values_intra=[1.1],
                 freq_cut_off=150, guess_found=False, mem='2GB', proc=2, max_cycles=30):
        """
        Initialize a TSOptimizer instance.

        Parameters:
        - rxn_id (int): Reaction ID.
        - reaction_smiles (str): SMILES representation of the reaction.
        - xtb_external_path (str): Path to external xtb executable.
        - xtb_solvent (str, optional): Solvent information for xTB calculations.
        - dft_solvent (str, optional): Solvent information for DFT calculations
        - reactive_complex_factor_values_inter (list, optional): Reactive complex
          factor values for intermolecular interactions.
        - reactive_complex_factor_values_intra (list, optional): Reactive complex
          factor values for intramolecular interactions.
        - freq_cut_off (int, optional): Frequency cutoff value.
        - guess_found (bool, optional): Indicates if a guess is found.
        - mem (str, optional): Gaussian memory specification. Defaults to '2GB'.
        - proc (int, optional): Number of processors. Defaults to 2.
        - max_cycles (int, optional): Maximal number of cycles in TS geometry search. Defaults to 30.
        """
        self.rxn_id = rxn_id
        self.reactant_smiles = reaction_smiles.split('>>')[0]
        self.product_smiles = reaction_smiles.split('>>')[-1]
        self.xtb_solvent = xtb_solvent
        self.dft_solvent = dft_solvent
        self.reactive_complex_factor_values_inter = reactive_complex_factor_values_inter
        self.reactive_complex_factor_values_intra = reactive_complex_factor_values_intra
        self.freq_cut_off = freq_cut_off
        self.xtb_external_path = xtb_external_path
        self.mem = mem
        self.proc = proc
        self.max_cycles = max_cycles

        self.charge, self.multiplicity = self.get_charge_and_multiplicity()

        self.reaction_dir = self.make_work_dir()

        self.g16_dir = self.make_sub_dir(sub_dir_name='g16_dir')
        self.rp_geometries_dir = self.make_sub_dir(sub_dir_name='rp_geometries')
        if not guess_found:
            self.path_dir = self.make_sub_dir(sub_dir_name='path_dir')
            self.ts_guesses_dir = self.make_sub_dir(sub_dir_name='preliminary_ts_guesses')
            self.final_guess_dir = self.make_sub_dir(sub_dir_name='final_ts_guess')
        else:
            self.path_dir = None
            self.final_guess_dir = None

        self.ts_guess_list = None

    def determine_ts(self, xtb=True, method='UB3LYP', basis_set='6-31G**'):
        """
        Determine the transition state for the given optimization parameters.

        Parameters:
        - xtb (bool, optional): Use external XTB for optimization.
        - method (str, optional): Level of theory for DFT optimization.
        - basis_set (str, optional): Basis set for DFT optimization.

        Returns:
        - bool: True if transition state is found, False otherwise.
        """
        if self.ts_guess_list is None:
            return False

        if xtb:
            if self.xtb_solvent is not None:
                method = f'external="{self.xtb_external_path} alpb={self.xtb_solvent}"'
            else:
                method = f'external="{self.xtb_external_path}"'
            basis_set = ''

        if self.dft_solvent is not None and not xtb:
            extra_commands = f'opt=(ts, calcall, noeigen, nomicro, MaxCycles={self.max_cycles}) SCRF=(Solvent={self.dft_solvent}, smd)'
        else:
            extra_commands = f'opt=(ts, calcall, noeigen, nomicro, MaxCycles={self.max_cycles})'

        remove_files_in_directory(self.g16_dir)

        for i, guess_file in enumerate(self.ts_guess_list):
            ts_search_inp_file = self.generate_g16_input_ts_opt(
                i, guess_file, method=method, basis_set=basis_set, extra_commands=extra_commands)
            os.chdir(self.reaction_dir)
            log_file = run_g16_ts_optimization(ts_search_inp_file)
            success = self.confirm_opt_transition_state(log_file, xtb, method, basis_set)

            if success:
                xyz_file = f'{os.path.splitext(log_file)[0]}.xyz'
                self.save_final_ts_guess_files(xyz_file, log_file)
                return True

        return False

    def set_ts_guess_list(self, reactive_complex_factor):
        """
        Set the list of transition state guesses based on the given reactive complex factor.

        Parameters:
        - reactive_complex_factor: Reactive complex factor value.
        """
        path = self.set_up_path_generator(reactive_complex_factor)
        ts_guess_list = self.obtain_ts_guesses_for_given_reactive_complex_factor(path)[:5]

        if ts_guess_list is not None:
            self.save_ts_guesses(ts_guess_list)

        self.ts_guess_list = ts_guess_list

    def modify_ts_guess_list(self, xyz_file_list):
        """
        Modify the list of transition state guess list if you already have a good guess.

        Parameters:
        - xyz_file_list (list): List of XYZ file paths for transition state guesses.
        """
        self.ts_guess_list = xyz_file_list


    def set_up_path_generator(self, reactive_complex_factor, n_conf=100):
        """
        Set up a path generator for the given reactive complex factor.

        Parameters:
        - reactive_complex_factor: Reactive complex factor value.
        - n_conf (int, optional): Number of conformations.

        Returns:
        - PathGenerator: Instance of the path generator.
        """
        # if intermolecular, immediately set up the PathGenerator object
        if '.' in self.reactant_smiles:
            path = PathGenerator(
                self.reactant_smiles, self.product_smiles, self.rxn_id, self.path_dir, self.rp_geometries_dir,
                self.xtb_solvent, reactive_complex_factor, self.freq_cut_off, self.charge, self.multiplicity, n_conf=n_conf
            )
        else:
            # Quick run to see if inversion is needed
            path = PathGenerator(
            self.reactant_smiles, self.product_smiles, self.rxn_id, self.path_dir, self.rp_geometries_dir,
            self.xtb_solvent, reactive_complex_factor, self.freq_cut_off, self.charge, self.multiplicity, n_conf=1
            )
            if len(path.formed_bonds) < len(path.broken_bonds) and '.' not in self.reactant_smiles:
                path = PathGenerator(
                    self.product_smiles, self.reactant_smiles, self.rxn_id, self.path_dir, self.rp_geometries_dir,
                    self.xtb_solvent, reactive_complex_factor, self.freq_cut_off, self.charge, self.multiplicity, n_conf=n_conf
                )
            else:
                path = PathGenerator(
                    self.reactant_smiles, self.product_smiles, self.rxn_id, self.path_dir, self.rp_geometries_dir,
                    self.xtb_solvent, reactive_complex_factor, self.freq_cut_off, self.charge, self.multiplicity, n_conf=n_conf
                )

        return path

    def obtain_ts_guesses_for_given_reactive_complex_factor(self, path):
        """
        Obtain transition state guesses for the given reactive complex factor.

        Parameters:
        - path: Path generator instance.

        Returns:
        - list or None: List of transition state guess files if found, None otherwise.
        """
        for _ in range(5):
            energies, potentials, path_xyz_files = path.get_path()
            if energies is not None:
                true_energies = list(np.array(energies) - np.array(potentials))
                guesses_list = self.determine_and_filter_local_maxima(true_energies, path_xyz_files, path.charge, path.multiplicity, path.solvent)
                if len(guesses_list) > 0:
                    return guesses_list

        return None

    def determine_and_filter_local_maxima(self, true_energies, path_xyz_files, charge, multiplicity, solvent):
        """
        Determine and filter local maxima in the path.

        Parameters:
        - true_energies (list): List of true energy values.
        - path_xyz_files (list): List of path XYZ files.
        - charge (int): Charge information.
        - multiplicity (int): Multiplicity information.

        Returns:
        - list: List of ranked transition state guess files based on energy.
        """
        # Find local maxima in path
        indices_local_maxima = find_local_max_indices(true_energies)

        # Validate the local maxima and store their energy values
        ts_guess_dict = {}
        idx_local_maxima_correct_mode = []
        for index in indices_local_maxima:
            ts_guess_file, _ = validate_ts_guess(path_xyz_files[index], self.reaction_dir, self.freq_cut_off, charge, multiplicity, solvent)
            if ts_guess_file is not None:
                idx_local_maxima_correct_mode.append(index)
                ts_guess_dict[ts_guess_file] = true_energies[index]

        # Sort guesses based on energy
        sorted_guess_dict = sorted(ts_guess_dict.items(), key=lambda x: x[1], reverse=True)
        ranked_guess_files = [item[0] for item in sorted_guess_dict]

        return ranked_guess_files

    def make_work_dir(self):
        """
        Create and return the main working directory.

        Returns:
        - str: Path to the main working directory.
        """
        dir_name = f'reaction_{self.rxn_id}'
        if dir_name in os.listdir(os.getcwd()):
            shutil.rmtree(os.path.join(os.getcwd(), dir_name))
        os.makedirs(os.path.join(os.getcwd(), dir_name))

        return os.path.join(os.getcwd(), dir_name)

    def make_sub_dir(self, sub_dir_name):
        """
        Create and return a subdirectory within the reaction directory.

        Parameters:
        - sub_dir_name (str): Name of the subdirectory.

        Returns:
        - str: Path to the created subdirectory.
        """
        if sub_dir_name in os.listdir(self.reaction_dir):
            shutil.rmtree(os.path.join(self.reaction_dir, sub_dir_name))
        os.makedirs(os.path.join(self.reaction_dir, sub_dir_name))

        return os.path.join(self.reaction_dir, sub_dir_name)

    def save_ts_guesses(self, ts_guesses_list):
        """
        Save transition state guess files to the ts_guesses_dir.

        Parameters:
        - ts_guesses_list (list): List of transition state guess files.
        """
        for i, ts_guess_file in enumerate(ts_guesses_list):
            shutil.copy(ts_guess_file, self.ts_guesses_dir)
            os.rename(
                os.path.join(self.ts_guesses_dir, os.path.basename(ts_guess_file)),
                os.path.join(self.ts_guesses_dir, f'ts_guess_{i}.xyz')
            )
    
    def generate_g16_input_ts_opt(self, idx, file_name, method, basis_set='',
                               extra_commands='opt=(ts, calcall, noeigen, nomicro, MaxCycles=30)'):
        """
        Generate Gaussian 16 input file for transition state optimization.

        Parameters:
        - idx (int): Index of the transition state guess.
        - file_name (str): Name of the transition state guess file.
        - method (str): Level of theory for optimization.
        - basis_set (str, optional): Basis set for optimization.
        - extra_commands (str, optional): Additional Gaussian commands.

        Returns:
        - str: Path to the generated Gaussian input file.
        """
        ts_search_inp_file = os.path.join(self.g16_dir, f'ts_guess_{idx}.com')
        xyz_to_gaussian_input(
            os.path.join(self.path_dir, file_name), ts_search_inp_file,
            method=method, basis_set=basis_set, extra_commands=extra_commands,
            charge=self.charge, multiplicity=self.multiplicity, mem=self.mem, 
            proc=self.proc 
        )

        return ts_search_inp_file

    def confirm_opt_transition_state(self, log_file, xtb=True, method='UB3LYP', basis_set='6-31G**'):
        """
        Confirm the optimized transition state by performing IRC calculations.

        Parameters:
        - log_file (str): Path to the Gaussian log file.
        - xtb (bool, optional): Use external XTB for IRC calculations.
        - method (str, optional): Level of theory for IRC calculations.
        - basis_set (str, optional): Basis set for IRC calculations.

        Returns:
        - bool: True if the optimized transition state is confirmed, False otherwise.
        """
        try:
            extract_transition_state_geometry(log_file, f'{os.path.splitext(log_file)[0]}.xyz')
            if xtb:
                irc_input_file_f, irc_input_file_r = generate_gaussian_irc_input(
                    f'{os.path.splitext(log_file)[0]}.xyz',
                    output_prefix=f'{os.path.splitext(log_file)[0]}_irc',
                    method=self.xtb_external_path,
                    mem=self.mem, proc=self.proc,
                    solvent=self.xtb_solvent, charge=self.charge, multiplicity=self.multiplicity
                )
            else:
                irc_input_file_f, irc_input_file_r = generate_gaussian_irc_input(
                    f'{os.path.splitext(log_file)[0]}.xyz',
                    output_prefix=f'{os.path.splitext(log_file)[0]}_irc',
                    method=f'{method}/{basis_set} ',
                    mem=self.mem, proc=self.proc,
                    solvent=self.dft_solvent, charge=self.charge, multiplicity=self.multiplicity
                )

            run_irc(irc_input_file_f)
            run_irc(irc_input_file_r)
            extract_irc_geometries(f'{os.path.splitext(irc_input_file_f)[0]}.log',
                               f'{os.path.splitext(irc_input_file_r)[0]}.log')

            reaction_correct = compare_molecules_irc(
                f'{os.path.splitext(irc_input_file_f)[0]}.xyz',
                f'{os.path.splitext(irc_input_file_r)[0]}.xyz',
                os.path.join(self.rp_geometries_dir, 'reactants_geometry.xyz'),
                os.path.join(self.rp_geometries_dir, 'products_geometry.xyz'),
                self.charge, self.multiplicity, self.xtb_solvent
            )

            return reaction_correct

        except Exception:
            return False

    def reaction_is_intramolecular(self):
        """
        Check if the reaction is intramolecular based on the number of reactant components.

        Returns:
        - bool: True if the reaction is intramolecular, False otherwise.
        """
        path = self.set_up_path_generator(reactive_complex_factor=1.2, n_conf=1)
        if len(path.reactant_smiles.split('.')) == 1:
            return True
        else:
            return False

    # TODO: What about triplet states?
    def get_charge_and_multiplicity(self):
        """
        Get the charge and multiplicity of the reactant molecule.

        Returns:
        - tuple: (charge, multiplicity)
        """
        mol = Chem.MolFromSmiles(self.reactant_smiles, ps)
        charge = Chem.GetFormalCharge(mol)
        total_electrons = 0

        for atom in mol.GetAtoms():
            # Add the atomic number
            total_electrons += atom.GetAtomicNum()

        # subtract the net charge
        total_electrons -= charge
        multiplicity = total_electrons % 2 + 1

        return charge, multiplicity

    def save_final_ts_guess_files(self, xyz_file, log_file):
        """
        Save the final transition state guess files to the final_guess_dir.

        Parameters:
        - xyz_file (str): Path to the XYZ file.
        - log_file (str): Path to the log file.

        Returns:
        - None
        """
        shutil.copy(xyz_file, self.final_guess_dir)
        shutil.copy(log_file, self.final_guess_dir)


def find_local_max_indices(numbers):
    """
    Find indices of local maxima in a list of numbers.

    Parameters:
    - numbers (list): List of numbers.

    Returns:
    - list: List of indices corresponding to local maxima.
    """
    local_max_indices = []
    for i in range(len(numbers) - 2, 0, -1):
        if numbers[i] > numbers[i - 1] and numbers[i] > numbers[i + 1]:
            local_max_indices.append(i)
    return local_max_indices
