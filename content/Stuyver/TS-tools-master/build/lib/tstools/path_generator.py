from rdkit import Chem
import numpy as np
import autode as ade
import os
from autode.conformers import conf_gen
from autode.conformers import conf_gen, Conformer
from scipy.spatial import distance_matrix
import copy
import subprocess
import shutil
import random
from itertools import product

from typing import Optional
from autode.smiles.smiles import init_smiles

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from tstools.utils import write_xyz_file_from_ade_atoms, NotConverged

ps = Chem.SmilesParserParams()
ps.removeHs = False
bohr_ang = 0.52917721090380

xtb = ade.methods.XTB()

metal_list = ['Al', 'Sb', 'Ag', 'As', 'Ba', 'Be', 'Bi', 'Cd', 'Ca', 'Cr', 'Co', 'Cu', 'Au', 'Fe', 
              'Pb', 'Mg', 'Mn', 'Hg', 'Mo', 'Ni', 'Pd', 'Pt', 'K', 'Rh', 'Rb', 'Ru', 'Sc', 'Ag', 
              'Na', 'Sr', 'Ta', 'Tl', 'Th', 'Ti', 'U', 'V', 'Y', 'Zn', 'Zr']

class PathGenerator:
    # Constants
    FC_CRUDE_LOWER_BOUND = 0.1
    FC_CRUDE_UPPER_BOUND = 4.0
    FC_CRUDE_INCREMENT = 0.1
    FC_CRUDE_ATTEMTPS = 5

    FC_REFINED_LOWER_BOUND = 0.09
    FC_REFINED_UPPER_BOUND = 0.03
    FC_REFINED_INCREMENT = 0.01
    FC_REFINED_ATTEMPTS = 2

    MIN_FC_LOWER_BOUND = 0.009
    MIN_FC_UPPER_BOUND = 0.005
    FC_INCREMENT = 0.001
    MAX_FORCE_CONSTANT = 0.1
    POTENTIAL_THRESHOLD = 0.005

    STRETCH_FACTOR_LOWER_BOUND = 1.0
    STRETCH_FACTOR_UPPER_BOUND = 1.3

    def __init__(self, reactant_smiles, product_smiles, rxn_id, path_dir, rp_geometries_dir, 
                 solvent=None, reactive_complex_factor=2.0, freq_cut_off=150, charge=0, multiplicity=1, n_conf=100):
        """
        Initialize a PathGenerator object.

        Parameters:
        - reactant_smiles (str): SMILES representation of the reactant.
        - product_smiles (str): SMILES representation of the product.
        - rxn_id (str): Identifier for the reaction.
        - path_dir (str): Path to the directory for storing generated paths.
        - rp_geometries_dir (str): Path to the directory containing reactant and product geometries.
        - solvent (str, optional): Solvent information.
        - reactive_complex_factor (float, optional): Reactive complex factor.
        - freq_cut_off (int, optional): Frequency cutoff.
        - charge (int, optional): Molecular charge.
        - multiplicity (int, optional): Molecular multiplicity.
        - n_conf (int, optional): Number of conformers.

        Returns:
        None
        """
        self.reactant_smiles = reactant_smiles
        self.product_smiles = product_smiles
        self.rxn_id = rxn_id
        self.path_dir = path_dir
        self.rp_geometries_dir = rp_geometries_dir
        self.solvent = solvent
        self.reactive_complex_factor = reactive_complex_factor
        self.freq_cut_off = freq_cut_off
        self.charge = charge
        self.multiplicity = multiplicity
        self.n_conf = n_conf

        os.chdir(self.path_dir)

        self.reactant_rdkit_mol = Chem.MolFromSmiles(reactant_smiles, ps)
        self.product_rdkit_mol = Chem.MolFromSmiles(product_smiles, ps)
        self.formed_bonds, self.broken_bonds = self.get_active_bonds_from_mols()

        self.atom_map_dict = {atom.GetAtomMapNum(): atom.GetIdx() for atom in self.reactant_rdkit_mol.GetAtoms()}
        self.atom_idx_dict = {atom.GetIdx(): atom.GetAtomMapNum() for atom in self.reactant_rdkit_mol.GetAtoms()}

        self.owning_dict_rsmiles = get_owning_mol_dict(reactant_smiles)
        self.owning_dict_psmiles = get_owning_mol_dict(product_smiles)

        self.reaction_is_organometallic = self.check_if_reaction_organometallic()

        self.formation_constraints = self.get_optimal_distances()

        self.stereo_correct_conformer_name = self.get_stereo_correct_conformer_name(n_conf)

        if self.n_conf > 1:
            self.minimal_fc = self.determine_minimal_fc()

    def get_path(self):
        """
        Attempt to generate a reaction path.

        Returns:
        Tuple or None: A tuple containing energies, potentials, and path XYZ files if successful,
                       or None if no valid path is found.
        """
        path_xyz_files = None

        if self.minimal_fc is not None:
            i = 0
            # Overstretch range a bit because otherwise, you may end up aborting the search prematurely
            for fc in np.arange(self.minimal_fc - PathGenerator.MIN_FC_LOWER_BOUND,
                                self.minimal_fc + PathGenerator.MIN_FC_UPPER_BOUND,
                                PathGenerator.FC_INCREMENT):
                
                # Only reinitialize the reactive complex after a first attempt
                if i == 0:
                    reactive_complex_xyz_file = f'{self.stereo_correct_conformer_name}_opt.xyz'
                else:
                    reactive_complex_xyz_file = self.get_reactive_complex(min(fc, PathGenerator.MAX_FORCE_CONSTANT))

                # Perform the biased optimization
                energies, coords, atoms, potentials = self.get_path_for_biased_optimization(reactive_complex_xyz_file, fc)

                # Check conditions for a valid path
                if potentials[-1] > min(fc, PathGenerator.POTENTIAL_THRESHOLD):
                    continue  # Haven't reached the products at the end of the biased optimization
                else:
                    if not (self.endpoint_is_product(atoms, coords) and self.beginpoint_is_reactant(atoms, coords)):
                        # If third time incorrect begin-/endpoint is reached, abort
                        if i > 2:
                            print(f'Incorrect reactant/product formed for {self.rxn_id}')
                            path_xyz_files = get_path_xyz_files(atoms, coords, fc)
                            return None, None, None
                        else:
                            i += 1
                            continue
                    else:
                        path_xyz_files = get_path_xyz_files(atoms, coords, fc)
                        self.save_rp_geometries(atoms, coords)
                        return energies, potentials, path_xyz_files

        return None, None, None

    def determine_minimal_fc(self):
        """
        Determine the minimal force constant for the reaction path.

        Returns:
        float or None: The refined minimal force constant, or None if not found.
        """
        minimal_fc_crude = self.screen_fc_range(
            PathGenerator.FC_CRUDE_LOWER_BOUND,
            PathGenerator.FC_CRUDE_UPPER_BOUND,
            PathGenerator.FC_CRUDE_INCREMENT,
            n_attempts=PathGenerator.FC_CRUDE_ATTEMTPS
        )

        if minimal_fc_crude is not None:
            # Overstretch range a bit because otherwise you may end up aborting the search prematurely
            minimal_fc_refined = self.screen_fc_range(
                minimal_fc_crude - PathGenerator.FC_REFINED_LOWER_BOUND,
                minimal_fc_crude + PathGenerator.FC_REFINED_UPPER_BOUND,
                PathGenerator.FC_REFINED_INCREMENT,
                n_attempts= PathGenerator.FC_REFINED_ATTEMPTS  # Only one retry for the refined range
            )

            return minimal_fc_refined

        return None
    
    def screen_fc_range(self, start, end, interval, n_attempts):
        """
        Screen a force constant range to find a suitable force constant.

        Parameters:
        - start (float): Starting value of the force constant range.
        - end (float): Ending value of the force constant range.
        - interval (float): Interval between force constant values.
        - n_attempts (int): Number of attempts to find a suitable force constant.

        Returns:
        float or None: The found force constant or None if not found.
        """
        if n_attempts > 2:
            self.stereo_correct_conformer_name = self.get_stereo_correct_conformer_name(self.n_conf)

        for fc in np.arange(start, end, interval):
            for _ in range(n_attempts):
                reactive_complex_xyz_file = self.get_reactive_complex(min(fc, PathGenerator.MAX_FORCE_CONSTANT))
                _, _, _, potentials = self.get_path_for_biased_optimization(reactive_complex_xyz_file, fc)

                if potentials[-1] < PathGenerator.POTENTIAL_THRESHOLD:
                    return fc
                else:
                    continue

        return None

    def get_formation_constraints_stretched(self):
        """
        Get stretched formation constraints for bonds that are to be stretched.

        Returns:
        dict: Dictionary of stretched formation constraints.
        """
        formation_constraints_to_stretch = self.get_bonds_to_stretch()
        formation_constraints_stretched = {}

        for bond, original_distance in self.formation_constraints.items():
            if bond in formation_constraints_to_stretch:
                stretch_factor = random.uniform(
                    PathGenerator.STRETCH_FACTOR_LOWER_BOUND * self.reactive_complex_factor,
                    PathGenerator.STRETCH_FACTOR_UPPER_BOUND * self.reactive_complex_factor
                )
                formation_constraints_stretched[bond] = stretch_factor * original_distance

        return formation_constraints_stretched
    
    def get_stereo_correct_conformer_name(self, n_conf=100):
        """
        Get the name of a stereochemistry-correct conformer.

        Parameters:
        - n_conf (int): Number of conformers to generate.

        Returns:
        str or None: The name of the conformer or None if not found.
        """
        if self.reactive_complex_factor > 0.01:
            formation_constraints_stretched = self.get_formation_constraints_stretched()
        else:
            formation_constraints_stretched = {}

        get_conformer_with_ade(self.reactant_smiles, self.reactant_rdkit_mol)
    
        stereochemistry_smiles = find_stereocenters(self.reactant_rdkit_mol)
        # only consider defined smiles stereocenters!
        stereo_elements_to_consider = {str(d['position']) for d in stereochemistry_smiles}
        # convert stereo_elements in strings, so that you can later on convert list into set
        stereochemistry_smiles = [str(d) for d in stereochemistry_smiles]

        ade_mol = ade.Molecule(f'input_reactants.xyz', charge=self.charge, mult=self.multiplicity)

        for node in ade_mol.graph.nodes:
            ade_mol.graph.nodes[node]['stereo'] = False

        bonds = []
        for bond in self.reactant_rdkit_mol.GetBonds():
            i,j = self.atom_map_dict[bond.GetBeginAtom().GetAtomMapNum()], self.atom_map_dict[bond.GetEndAtom().GetAtomMapNum()]

            if (i, j) not in formation_constraints_stretched and (j, i) not in formation_constraints_stretched:
                bonds.append((i,j))

        ade_mol.graph.edges = bonds

        # find good starting conformer
        for n in range(n_conf):
            atoms = conf_gen.get_simanl_atoms(species=ade_mol, dist_consts=formation_constraints_stretched, conf_n=n, save_xyz=False) # set save_xyz to false to ensure new optimization
            conformer = Conformer(name=f"conformer_reactants_init", atoms=atoms, charge=self.charge, dist_consts=formation_constraints_stretched)
            write_xyz_file_from_ade_atoms(atoms, f'{conformer.name}.xyz')
            stereochemistry_conformer = get_stereochemistry_from_conformer_xyz(f'{conformer.name}.xyz', self.reactant_smiles)
            stereochemistry_conformer = [str(d) for d in stereochemistry_conformer if str(d['position']) in stereo_elements_to_consider]

            if set(stereochemistry_smiles) == set(stereochemistry_conformer):
                return conformer.name

        # print that there is an error with the stereochemistry only when you do a full search, i.e., n_conf > 1
        if n_conf > 1:
            print(f'No stereo-compatible conformer found for reaction {self.rxn_id}')

        return conformer.name

    def get_reactive_complex(self, fc):
        """
        Get the optimized reactive complex.

        Parameters:
        - force_constant (float): The force constant for the optimization.

        Returns:
        str: Path to the optimized reactive complex XYZ file.
        """
        if self.reactive_complex_factor > 0.01:
            formation_constraints_stretched = self.get_formation_constraints_stretched()
        else:
            formation_constraints_stretched = {}

        self.optimize_reactive_complex(formation_constraints_stretched, fc)

        return f'{self.stereo_correct_conformer_name}_opt.xyz'

    def optimize_reactive_complex(self, formation_constraints_stretched, fc):
        """
        Optimize the geometry of a reactive complex using the xTB quantum chemistry package.

        Args:
            formation_constraints_stretched (dict): A dictionary containing the formation constraints stretched.
            fc (float): Force constant for the constraint during optimization.

        Raises:
            RuntimeError: If an error occurs during the xTB optimization process.
        """
        xtb_input_path = f'{self.stereo_correct_conformer_name}.inp'

        with open(xtb_input_path, 'w') as f:
            f.write('$constrain\n')
            f.write(f'    force constant={fc}\n')
            for key, val in formation_constraints_stretched.items():
                f.write(f'    distance: {key[0] + 1}, {key[1] + 1}, {val}\n')
            f.write('$end\n')
        cmd = f'xtb {self.stereo_correct_conformer_name}.xyz --opt --input {xtb_input_path} -v --charge {self.charge} '

        if self.solvent is not None:
            cmd += f'--alpb {self.solvent} '
        if self.multiplicity != 1:
            cmd += f'--uhf {self.multiplicity - 1} '

        try:
            with open(f'{self.stereo_correct_conformer_name}.out', 'w') as out:
                process = subprocess.Popen(cmd.split(), stderr=subprocess.DEVNULL, stdout=out)
                process.wait()
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f'Error during XTB optimization: {e}')
        
        try:
            os.rename('xtbopt.xyz', f'{self.stereo_correct_conformer_name}_opt.xyz')
        except: # xTB calculation didn't converge
            raise NotConverged(self.rxn_id)

    def get_path_for_biased_optimization(self, reactive_complex_xyz_file, force_constant):
        """
        Perform biased optimization and retrieve valid path information.

        Parameters:
        - reactive_complex_xyz_file (str): Path to the reactive complex XYZ file.
        - force_constant (float): The force constant for the optimization.

        Returns:
        Tuple: Tuple containing valid energies, coordinates, atoms, and potentials.
        """
        log_file = self.xtb_optimize_with_applied_potentials(reactive_complex_xyz_file, force_constant)
        all_energies, all_coords, all_atoms = read_energy_coords_file(log_file)

        valid_energies, valid_coords, valid_atoms = [], [], []
        for i, coords in enumerate(all_coords):
            valid_coords.append(coords)
            valid_atoms.append(all_atoms[i])
            valid_energies.append(all_energies[i])

        potentials = determine_potential(valid_coords, self.formation_constraints, force_constant)

        return valid_energies, valid_coords, valid_atoms, potentials

    def xtb_optimize_with_applied_potentials(self, reactive_complex_xyz_file, fc):
        """
        Perform XTB optimization with applied potentials.

        Parameters:
        - reactive_complex_xyz_file (str): Path to the reactive complex XYZ file.
        - force_constant (float): The force constant for the applied potentials.

        Returns:
        str: Path to the XTB optimization log file.
        """
        xtb_input_path = f'{os.path.splitext(reactive_complex_xyz_file)[0]}.inp'

        with open(xtb_input_path, 'w') as f:
            f.write('$constrain\n')
            f.write(f'    force constant={fc}\n')
            for key, val in self.formation_constraints.items():
                f.write(f'    distance: {key[0] + 1}, {key[1] + 1}, {val}\n')
            f.write('$end\n')

        cmd = f'xtb {reactive_complex_xyz_file} --opt --input {xtb_input_path} -v --charge {self.charge} '

        if self.solvent is not None:
            cmd += f'--alpb {self.solvent} '
        if self.multiplicity != 1:
            cmd += f'--uhf {self.multiplicity - 1} '

        try:
            with open(f'{os.path.splitext(reactive_complex_xyz_file)[0]}_path.out', 'w') as out:
                process = subprocess.Popen(cmd.split(), stderr=subprocess.DEVNULL, stdout=out)
                process.wait()
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f'Error during XTB optimization: {e}')

        try:
            os.rename('xtbopt.log', f'{os.path.splitext(reactive_complex_xyz_file)[0]}_path.log')
        except: # xTB calculation didn't converge
            raise NotConverged(self.rxn_id)

        return f'{os.path.splitext(reactive_complex_xyz_file)[0]}_path.log'     

    def get_bonds_to_stretch(self):
        """
        Get the set of bonds to be stretched based on formation constraints.

        Returns:
        set: Set of bonds to be stretched.
        """
        bonds_to_stretch = set()

        for bond in self.formation_constraints.keys():
            atom_1, atom_2 = bond
            owner_1 = self.owning_dict_rsmiles[self.atom_idx_dict[atom_1]]
            owner_2 = self.owning_dict_rsmiles[self.atom_idx_dict[atom_2]]

            if owner_1 != owner_2:
                bonds_to_stretch.add(bond)

        # If no intermolecular bonds found, consider all bonds
        if not bonds_to_stretch:
            bonds_to_stretch = set(self.formation_constraints.keys())

        return bonds_to_stretch

    def get_active_bonds_from_mols(self):
        """
        Identify formed and broken bonds between reactant and product molecules.

        Returns:
        set: Formed bonds in the product.
        set: Broken bonds in the reactant.
        """
        reactant_bonds = get_bonds(self.reactant_rdkit_mol)
        product_bonds = get_bonds(self.product_rdkit_mol)

        formed_bonds = product_bonds - reactant_bonds
        broken_bonds = reactant_bonds - product_bonds

        return formed_bonds, broken_bonds
    
    def get_optimal_distances(self):
        """
        Calculate optimal distances for formed bonds in the product 
        (add additional distance constraints if organometallic system).

        Returns:
        dict: Dictionary of optimal distances for formed bonds.
        """
        optimal_distances = {}
        product_smiles = [smi for smi in self.product_smiles.split('.')]
        product_molecules = [Chem.MolFromSmiles(smi, ps) for smi in self.product_smiles.split('.')]
        formed_bonds = self.formed_bonds

        atoms_involved_in_formed_bonds = []

        for bond in formed_bonds:
            atom_i = int(bond[0])
            atom_j = int(bond[1])

            idx1, idx2 = self.atom_map_dict[atom_i], self.atom_map_dict[atom_j]

            mol, mol_dict, smiles = self.get_mol_and_mol_dict(atom_i, atom_j, product_molecules, product_smiles)
            current_bond_length = self.obtain_current_distance(mol, mol_dict, smiles, atom_i, atom_j)
 
            optimal_distances[idx1, idx2] = current_bond_length

            # for metal-containing bonds, add the atoms that involve main group elements to the atoms_involved_in_formed_bonds list
            if mol.GetAtomWithIdx(mol_dict[atom_i]).GetSymbol() not in metal_list and \
                  mol.GetAtomWithIdx(mol_dict[atom_j]).GetSymbol() in metal_list:
                atoms_involved_in_formed_bonds.append(atom_i)
            if mol.GetAtomWithIdx(mol_dict[atom_i]).GetSymbol() in metal_list and \
                  mol.GetAtomWithIdx(mol_dict[atom_j]).GetSymbol() not in metal_list:
                atoms_involved_in_formed_bonds.append(atom_j)

        if self.reaction_is_organometallic == True:
            for atom_i, atom_j in list(product(atoms_involved_in_formed_bonds, repeat=2)):
                if (min(atom_i, atom_j), max(atom_i, atom_j)) in self.broken_bonds:
                    idx1, idx2 = self.atom_map_dict[atom_i], self.atom_map_dict[atom_j]
                    mol, mol_dict, smiles = self.get_mol_and_mol_dict(atom_i, atom_j, product_molecules, product_smiles)
                    current_distance = self.obtain_current_distance(mol, mol_dict, smiles, atom_i, atom_j)
                    optimal_distances[min(idx1, idx2), max(idx1, idx2)] = current_distance
                    break
                else:
                    continue

        return optimal_distances
    
    def get_mol_and_mol_dict(self, atom_i, atom_j, product_molecules, product_smiles):
        """
        Retrieve the molecule containing atoms 'atom_i' and 'atom_j' and create a dictionary
        mapping atom map numbers to atom indices in this molecule.

        Parameters:
        - atom_i (int): Index of the first atom in the bond.
        - atom_j (int): Index of the second atom in the bond.
        - product_molecules (dict): Dictionary mapping atom indices to corresponding molecules.

        Returns:
        - mol (rdkit.Chem.Mol): A deep copy of the molecule containing 'atom_i' and 'atom_j'.
        - mol_dict (dict): A dictionary mapping atom map numbers to atom indices in 'mol'.

        Raises:
        - KeyError: If atoms 'atom_i' and 'atom_j' belong to different molecules.
        """
        if self.owning_dict_psmiles[atom_i] == self.owning_dict_psmiles[atom_j]:
                mol = copy.deepcopy(product_molecules[self.owning_dict_psmiles[atom_i]])
                smiles = product_smiles[self.owning_dict_psmiles[atom_i]]
        else:
            raise KeyError("Atoms in the bond belong to different molecules.")

        mol_dict = {atom.GetAtomMapNum(): atom.GetIdx() for atom in mol.GetAtoms()}
        
        return mol, mol_dict, smiles
    
    def obtain_current_distance(self, mol, mol_dict, smiles, atom_i, atom_j):
        """
        Calculate the current distance between two atoms in the molecule.

        Parameters:
        - mol (rdkit.Chem.Mol): The molecule containing the atoms.
        - mol_dict (dict): A dictionary mapping atom map numbers to atom indices in 'mol'.
        - smiles (string): The molecule SMILES.
        - atom_i (int): Index of the first atom.
        - atom_j (int): Index of the second atom.

        Returns:
        - current_distance (float): The distance between the specified atoms in the molecule.

        Notes:
        - The function internally uses the obtain_dist_matrix method to compute the distance matrix.
        """  
        dist_matrix = self.obtain_dist_matrix(mol, smiles)
        current_distance = dist_matrix[mol_dict[atom_i], mol_dict[atom_j]]

        return current_distance

    def obtain_dist_matrix(self, mol, smiles):
        """
        Obtain the distance matrix for the atoms in the molecule.

        Parameters:
        - mol (rdkit.Chem.Mol): The molecule for which the distance matrix is to be calculated.
        - smiles (string): the molecule SMILES string.

        Returns:
        - dist_matrix (numpy.ndarray): The distance matrix representing pairwise distances between atoms.

        Notes:
        - The function modifies the atom map numbers in the molecule to ensure correct processing.
        - It temporarily writes a temporary XYZ file to avoid reordering of atoms by autodE.
        - The distance matrix is calculated using the autodE library.

        Raises:
        - Any exceptions raised during the optimization process using autodE.
        """
        [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
        get_conformer_with_ade(smiles, mol, output_file_name='tmp.xyz')
        charge = Chem.GetFormalCharge(mol)
        multiplicity = get_multiplicity(mol)

        if self.solvent is not None:
            ade_mol = ade.Molecule('tmp.xyz', name='tmp', charge=charge, mult=multiplicity, solvent_name=self.solvent)
        else:
            ade_mol = ade.Molecule('tmp.xyz', name='tmp', charge=charge, mult=multiplicity)

        ade_mol.conformers = [conf_gen.get_simanl_conformer(ade_mol)]
        ade_mol.conformers[0].optimise(method=xtb)
        dist_matrix = distance_matrix(ade_mol.conformers[0].coordinates, ade_mol.conformers[0].coordinates)
    
        return dist_matrix

    def save_rp_geometries(self, final_atoms, final_coords):
        """
        Save reactant and product geometries to XYZ files.

        Parameters:
        final_atoms (list): List of atoms for reactant and product.
        final_coords (list): List of coordinates for reactant and product.
        """
        reactant_xyz_path = os.path.join(self.rp_geometries_dir, 'reactants_geometry.xyz')
        product_xyz_path = os.path.join(self.rp_geometries_dir, 'products_geometry.xyz')

        write_xyz_file_from_atoms_and_coords(final_atoms[0], final_coords[0], reactant_xyz_path)
        write_xyz_file_from_atoms_and_coords(final_atoms[-1], final_coords[-1], product_xyz_path)

    def beginpoint_is_reactant(self, final_atoms, final_coords):
        """
        Check if the final geometry corresponds to the reactant.

        Parameters:
        final_atoms (list): List of atoms for the final geometry.
        final_coords (list): List of coordinates for the final geometry.

        Returns:
        bool: True if the final geometry corresponds to the reactant, False otherwise.
        """
        reactant_bonds = [(min(self.atom_map_dict[atom1], self.atom_map_dict[atom2]), max(self.atom_map_dict[atom1], self.atom_map_dict[atom2]))
                      for atom1, atom2 in get_bonds(self.reactant_rdkit_mol)]
        write_xyz_file_from_atoms_and_coords(final_atoms[0], final_coords[0], 'reactant_geometry.xyz')
        ade_mol_r = ade.Molecule('reactant_geometry.xyz', name='reactant_geometry', charge=self.charge, mult=self.multiplicity)

        # you need to treat organometallic reactions differently because non-bonded atoms may be close 
        # enough for distance-based bond assignment to be triggered
        if not self.reaction_is_organometallic:
            return set(ade_mol_r.graph.edges) == set(reactant_bonds)
        else:
            return set(ade_mol_r.graph.edges).issubset(set(reactant_bonds))
        
    def endpoint_is_product(self, final_atoms, final_coords):
        """
        Check if the final geometry corresponds to the product.

        Parameters:
        final_atoms (list): List of atoms for the final geometry.
        final_coords (list): List of coordinates for the final geometry.

        Returns:
        bool: True if the final geometry corresponds to the product, False otherwise.
        """
        product_bonds = [(min(self.atom_map_dict[atom1], self.atom_map_dict[atom2]), max(self.atom_map_dict[atom1], self.atom_map_dict[atom2]))
                     for atom1, atom2 in get_bonds(self.product_rdkit_mol)]
        write_xyz_file_from_atoms_and_coords(final_atoms[-1], final_coords[-1], 'products_geometry.xyz')
        ade_mol_p = ade.Molecule('products_geometry.xyz', name='products_geometry', charge=self.charge, mult=self.multiplicity)

        # you cannot do this check for organometallic reactions because non-bonded atoms may be close 
        # enough for distance based-bond assignment to be triggered
        if not self.reaction_is_organometallic:
            return set(ade_mol_p.graph.edges) == set(product_bonds)
        else:
            return set(ade_mol_p.graph.edges).issubset(set(product_bonds))
    
    def check_if_reaction_organometallic(self):
        """
        Check if the reactant molecule contains any organometallic atoms.

        Returns:
        - bool: True if organometallic atoms are present, False otherwise.

        Notes:
        - The function examines the atomic symbols of atoms in the reactant molecule.
        - It checks if any of the symbols match those in the 'metal_list'.
        """
        symbol_list = [atom.GetSymbol() for atom in self.reactant_rdkit_mol.GetAtoms()]

        for symbol in symbol_list:
            if symbol in metal_list:
                return True
            else:
                continue
        
        return False
    
def get_multiplicity(mol):
    """
    Get the multiplicity of a molecule.

    Returns:
    - int: multiplicity
    """
    charge = Chem.GetFormalCharge(mol)
    total_electrons = 0

    for atom in mol.GetAtoms():
        # Add the atomic number
        total_electrons += atom.GetAtomicNum()

    # subtract the net charge
    total_electrons -= charge
    multiplicity = total_electrons % 2 + 1

    return multiplicity


def get_owning_mol_dict(smiles):
    """
    Create a dictionary mapping atom map numbers to the index of the molecule to which they belong.

    Parameters:
    reaction_smiles (str): Reaction SMILES string.

    Returns:
    dict: A dictionary where keys are atom map numbers and values are the corresponding molecule indices.
    """
    molecules = [Chem.MolFromSmiles(smi, ps) for smi in smiles.split('.')]
    owning_mol_dict = {}

    for mol_index, mol in enumerate(molecules):
        for atom in mol.GetAtoms():
            owning_mol_dict[atom.GetAtomMapNum()] = mol_index

    return owning_mol_dict


def get_bonds(mol):
    """
    Get the bond strings of a molecule.

    Args:
        mol (Chem.Mol): Molecule.

    Returns:
        set: Set of bond strings.
    """
    bonds = set()
    for bond in mol.GetBonds():
        atom_1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetAtomMapNum()
        atom_2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomMapNum()

        if atom_1 < atom_2:
            bonds.add((atom_1, atom_2))
        else:
            bonds.add((atom_2, atom_1))

    return bonds


def read_energy_coords_file(file_path):
    """
    Read energy and coordinate information from a file.

    Args:
        file_path (str): The path to the file.

    Returns:
        Tuple: A tuple containing the energy values, coordinates, and atom symbols.
    """
    all_energies = []
    all_coords = []
    all_atoms = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            # read energy value from line starting with "energy:"
            if len(lines[i].split()) == 1 and lines[i+1].strip().startswith("energy:"):
                energy_line = lines[i+1].strip()
                energy_value = float(energy_line.split()[1])
                all_energies.append(energy_value)
                i += 2
            else:
                raise ValueError(f"Unexpected line while reading energy value: {energy_line}")
            # read coordinates and symbols for next geometry
            coords = []
            atoms = []
            while i < len(lines) and len(lines[i].split()) != 1:
                atoms.append(lines[i].split()[0])
                coords.append(np.array(list(map(float,lines[i].split()[1:]))))
                i += 1

            all_coords.append(np.array(coords))
            all_atoms.append(atoms)

    return np.array(all_energies), all_coords, all_atoms


def determine_potential(all_coords, constraints, force_constant):
    """
    Determine the potential energy for a set of coordinates based on distance constraints and a force constant.

    Args:
        all_coords (list): A list of coordinate arrays.
        constraints (dict): A dictionary specifying the atom index pairs and their corresponding distances.
        force_constant (float): The force constant to apply to the constraints.

    Returns:
        list: A list of potential energy values.
    """
    potentials = []
    for coords in all_coords:
        potential = 0
        dist_matrix = distance_matrix(coords, coords)
        for key, val in constraints.items():
            actual_distance = dist_matrix[key[0], key[1]] - val
            potential += force_constant * angstrom_to_bohr(actual_distance) ** 2
        potentials.append(potential)

    return potentials


def angstrom_to_bohr(distance_angstrom):
    """
    Convert distance in angstrom to bohr.

    Args:
        distance_angstrom (float): Distance in angstrom.

    Returns:
        float: Distance in bohr.
    """
    return distance_angstrom * 1.88973


def get_path_xyz_files(atoms, coords, force_constant):
    """
    Save a series of XYZ files representing the path along a reaction coordinate.

    Parameters:
    atoms (list): List of atom objects for each step in the path.
    coords (list): List of coordinate arrays for each step in the path.
    force_constant (float): The force constant applied to the path.

    Returns:
    list: List of filenames for the saved XYZ files.
    """
    path_xyz_files = []
    folder_name = f'path_xyzs_{force_constant:.4f}'
    if folder_name in os.listdir():
        shutil.rmtree(folder_name)
    os.makedirs(folder_name)

    for i in range(len(atoms)):
        filename = write_xyz_file_from_atoms_and_coords(
            atoms[i],
            coords[i],
                f'{folder_name}/path_{force_constant:.4}_{i}.xyz'
            )
        path_xyz_files.append(filename)  

    return path_xyz_files


def write_xyz_file_from_atoms_and_coords(atoms, coords, filename):
    """
    Write an XYZ file from a list of atoms and coordinates.

    Args:
        atoms: The list of atom symbols.
        coords: The list of atomic coordinates.
        filename: The name of the XYZ file to write.

    Returns:
        str: The name of the written XYZ file.
    """
    with open(filename, 'w') as f:
        f.write(f'{len(atoms)}\n')
        f.write("test \n")
        for i, coord in enumerate(coords):
            x, y, z = coord
            f.write(f"{atoms[i]} {x:.6f} {y:.6f} {z:.6f}\n")
    return filename


def find_stereocenters(mol):
    """
    Identify stereocenters (chirality and cis/trans bonds) in a molecule.

    Parameters:
    mol (Chem.Mol): RDKit molecule object.

    Returns:
    list: List of dictionaries representing stereocenters.
        Each dictionary contains keys:
        - 'type': Either 'chirality' or 'cis/trans'.
        - 'position': Atom map numbers or bond indices associated with the stereocenter.
        - 'descriptor': Chirality descriptor or cis/trans bond stereochemistry.
    """
    stereocenters = []

    # Find chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    for atom_idx, center_info in chiral_centers:
        stereocenter = {
            'type': 'chirality',
            'position': mol.GetAtomWithIdx(atom_idx).GetAtomMapNum(),
            'descriptor': center_info
        }
        stereocenters.append(stereocenter)

    # Find cis/trans bonds
    for bond in mol.GetBonds():
        if bond.GetStereo() > 0:
            stereocenter = {
                'type': 'cis/trans',
                'position': {bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()},
                'descriptor': bond.GetStereo()
            }
            stereocenters.append(stereocenter)

    return stereocenters


def get_stereochemistry_from_conformer_xyz(xyz_file, smiles):
    """
    Get stereochemistry information from an XYZ file.

    Args:
        xyz_file: The XYZ file.
        smiles: The SMILES string.

    Returns:
        object: The molecule with stereochemistry.
        list: The stereochemistry information.
    """
    mol = Chem.MolFromSmiles(smiles, ps)
    Chem.RemoveStereochemistry(mol)
    no_stereo_smiles = Chem.MolToSmiles(mol)
    mol = add_xyz_conformer(no_stereo_smiles, xyz_file)

    mol.GetConformer()
    Chem.AssignStereochemistryFrom3D(mol)
    stereochemistry = find_stereocenters(mol)

    return stereochemistry


def add_xyz_conformer(smiles, xyz_file):
    """
    Add an XYZ conformer to the molecule.

    Args:
        smiles: The SMILES string.
        xyz_file: The XYZ file.

    Returns:
        object: The molecule with the added conformer.
    """
    mol = Chem.MolFromSmiles(smiles, ps)
    
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
        num_atoms = int(lines[0])
        coords = []
        symbols = []
        for i in range(2, num_atoms+2):
            line = lines[i].split()
            symbol = line[0]
            x, y, z = map(float, line[1:])
            symbols.append(symbol)
            coords.append((x, y, z))

    conformer = Chem.Conformer(num_atoms)
    for i, coord in enumerate(coords):
        conformer.SetAtomPosition(i, coord)
    mol.AddConformer(conformer)
    
    return mol


# TODO: Are you sure that the ordering will now always be correct -> e.g., if you have two molecules, it could be that their order gets reversed in the SMILES, no? 
# Check this, and if not, try to assert that the ordering is the same
def get_conformer_with_ade(smiles, mol, output_file_name='input_reactants.xyz'):
    '''
    Generate conformers for a molecule specified by SMILES using autodE and save them in an XYZ file.

    Parameters:
        smiles (str): SMILES representation of the molecule.
        mol (ade.Molecule): ADE Molecule object representing the full molecule.
        output_file_name (str, optional): Name of the output XYZ file. Default is 'input_reactants.xyz'.

    Notes:
        - The function generates conformers for each component in the given SMILES string.
        - It creates temporary XYZ files for each conformer and combines them into a single XYZ file.
        - The combined XYZ file contains all conformers of the molecule along with the full molecule's structure.
    '''
    ind_mol_xyz_list = []

    for i, smi in enumerate(smiles.split('.')):
        ade_tmp_mol = ModifiedMolecule(name=f'tmp_{i}', smiles=smi)
        # just to ensure that autodE will later on generate a reasonable geometry from scratch
        if '.' in smiles:
            for atom in ade_tmp_mol.atoms:
                atom.coord = np.array([0.0, 0.0, 0.0])
        write_xyz_file_from_ade_atoms(ade_tmp_mol.atoms, f'tmp_{i}.xyz')
        ind_mol_xyz_list.append(f'tmp_{i}.xyz')

    combine_xyz_files(output_file_name, ind_mol_xyz_list, mol)


# you want to try custom initialization of the SMILES to ensure that you have the correct atom ordering
class ModifiedMolecule(ade.Molecule):
    '''
    A modified version of the molecule class from autodE where you always parse according to the ordering of the SMILES
    '''
    def _init_smiles(self, smiles: str, charge: Optional[int]):
        init_smiles(self, smiles)

        if charge is not None and charge != self._charge:
            raise ValueError(
                "SMILES charge was not the same as the "
                f"defined value. {self._charge} ≠ {charge}"
            )
        
        return None
    

def combine_xyz_files(output_file, input_files, full_mol):
    '''
    Merge contents of multiple XYZ files into a single output file.

    Parameters:
        input_files (list of str): List of paths to input XYZ files to be merged.
        output_file (str): Path to the output file where merged content will be written.
        full_mol (ade.Molecule): The full molecule object representing the combined system.

    Raises:
        FileNotFoundError: If any of the input files does not exist.
    
    Notes:
        - The function reads each input XYZ file, skipping the first two lines (assumed header).
        - It then appends the remaining lines from each input file to a list.
        - Finally, it writes the combined atomic count (number of atoms in full_mol) followed by the
          contents of the line_list to the output file.
    '''
    line_list = []

    with open(output_file, 'w') as output:
        for input_file in input_files:
            with open(input_file, 'r') as input:
                # Skip the first two lines (header)
                for _ in range(2):
                    next(input)
                # Copy the remaining lines to the output file
                for line in input:
                    line_list.append(line)

        output.write(f'{len(full_mol.GetAtoms())}\n\n')

        for line in line_list:
            output.write(line)


if __name__ == '__main__':
    shutil.rmtree('test')
    os.mkdir('test')
    os.chdir('test')
    #reactant_smiles = '[N+:1](=[B-:2](/[H:6])[H:7])(\[H:8])[H:9].[N+:3](=[B-:4](/[H:11])[H:12])(\[H:5])[H:10]'
    #product_smiles = '[N+:1]([B-:2]([H:6])([H:7])[H:12])([B:4]([N:3]([H:5])[H:10])[H:11])([H:8])[H:9]'
    #reactant_smiles = '[H:1]/[C:2](=[C:3](/[H:5])[O:6][H:7])[H:4].[O:8]=[C:9]([H:10])[H:11]'
    #product_smiles = '[H:1][C:2]([C:3]([H:5])=[O:6])([H:4])[C:9]([O:8][H:7])([H:10])[H:11]'
    reactant_smiles = '[H:1]/[C:2](=[C:3](/[H:5])[O:6][H:7])[H:4].[O:8]=[C:9]([H:10])[H:11]'
    product_smiles = '[H:1]/[C:2](=[C:3](\[O:6][H:7])[C:9]([O:8][H:5])([H:10])[H:11])[H:4]'
    reaction = PathGenerator(reactant_smiles, product_smiles, 'R1', 
                             '/Users/thijsstuyver/Desktop/reaction_profile_generator/path_test', 
                             '/Users/thijsstuyver/Desktop/reaction_profile_generator/rp_test')
    reaction.get_ts_guesses_from_path()
