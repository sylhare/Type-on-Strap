import numpy as np
import autode as ade
import os
from scipy.spatial import distance_matrix
import re
from autode.species.species import Species
from typing import Optional
import subprocess

xtb = ade.methods.XTB()


def validate_ts_guess(ts_guess_file, path, freq_cut_off=150, charge=0, multiplicity=1, solvent=None):
    """
    Validate the transition state (TS) guess file based on frequency and displacement information.

    Parameters:
    - ts_guess_file (str): Path to the TS guess file.
    - path (str): Path to the relevant data.
    - freq_cut_off (int, optional): Frequency cutoff for validation. Defaults to 150.
    - charge (int, optional): Charge of the system. Defaults to 0.
    - multiplicity (int, optional): Multiplicity of the system. Defaults to 1.
    - solvent (str, optional): Solvent information. Defaults to None.

    Returns:
    - Tuple[str, int] or Tuple[None, None]: If the TS guess is valid, returns the TS guess file
      and its frequency; otherwise, returns None for both.
    """
    # get all information about main imaginary mode
    freq, main_displacement_is_active = extract_info_ts_file(ts_guess_file, path, charge, multiplicity, solvent)

    if freq < -freq_cut_off and main_displacement_is_active:
        return ts_guess_file, freq
    else:
        return None, None


def extract_info_ts_file(ts_file, path, charge, multiplicity, solvent):
    """
    Extract information related to a transition state (TS) from a directory.

    Args:
        ts_file (str): The directory containing TS files.
        path (str): The path containing the reactant and product xyz-files.
        charge (int): The charge of the system.
        multiplicity (int): The multiplicity of the system.
        solvent (str): The name of the solvent.

    Returns:
        tuple: A tuple containing the following information:
            - float: Frequency of the TS.
            - dict: Active bonds involved in the imaginary mode, with bond indices as keys and displacement values as values.
            - dict: Extra bonds involved in the imaginary mode, with bond indices as keys and displacement values as values.
            - set: Active bonds forming during the TS.
            - set: Active bonds breaking during the TS.
            - numpy.ndarray: Distance matrix for reactant molecules.
            - numpy.ndarray: Distance matrix for product molecules.
            - numpy.ndarray: Distance matrix for the TS geometry.

    This function analyzes the provided TS directory to determine if it represents an imaginary mode and extracts various relevant information, 
    including bond displacements, active bonds forming and breaking, and distance matrices for the reactants, products, and TS geometry.

    Note:
    - The bond displacement cutoff (cut_off) is used to filter small bond displacements. 
        Bonds with displacements below this threshold are ignored.
    """
    # Obtain reactant, product, and transition state molecules
    reactant_file, product_file = get_xyzs(path)
    reactant, product, ts_mol = get_ade_molecules(reactant_file, product_file, ts_file, charge, multiplicity)   

    # Compute the displacement along the imaginary mode
    _ = get_negative_frequencies(ts_file, charge, solvent)
    normal_mode, freq = read_first_normal_mode('g98.out')
    f_displaced_species = displaced_species_along_mode(ts_mol, normal_mode, disp_factor=1)
    b_displaced_species = displaced_species_along_mode(reactant, normal_mode, disp_factor=-1)

    # Compute distance matrices -- TS geometry obtained through displacement along imaginary mode
    f_distances = distance_matrix(f_displaced_species.coordinates, f_displaced_species.coordinates)
    b_distances = distance_matrix(b_displaced_species.coordinates, b_displaced_species.coordinates)

    # Compute delta_mode
    delta_mode = f_distances - b_distances

    # Get all the bonds in both reactants and products
    all_bonds = set(product.graph.edges).union(set(reactant.graph.edges))

    # Identify active forming and breaking bonds
    active_bonds_forming = set(product.graph.edges).difference(set(reactant.graph.edges))
    active_bonds_breaking = set(reactant.graph.edges).difference(set(product.graph.edges))
    active_bonds = active_bonds_forming.union(active_bonds_breaking)

    # Check if main bond displacement in mode corresponds to active bond
    displacement_dict = {}
    for bond in all_bonds:
        displacement_dict[bond] = abs(delta_mode[bond[0],bond[1]])

    max_displacement_bond = max(displacement_dict, key=displacement_dict.get)

    if max_displacement_bond in active_bonds:
        main_displacement_is_active = True
    else:
        main_displacement_is_active = False
    
    return freq, main_displacement_is_active


def read_first_normal_mode(filename):
    """
    Read the first normal mode from the specified file.

    Args:
        filename (str): The name of the file to read.

    Returns:
        numpy.ndarray: Array representing the normal mode.
        float: Frequency value.
    """
    normal_mode = []
    pattern = r'\s+(\d+)\s+\d+\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)'

    # Open the file and read its contents
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Iterate over the lines and find the matching pattern
        for line in lines:
            # Check if the line contains a frequency
            if 'Frequencies' in line:
                # Extract the frequency value from the line
                frequency = float(line.split('--')[1].split()[0])

                # Iterate over the lines below the frequency line
                for sub_line in lines[lines.index(line) + 7:]:
                    # Check if the line matches the pattern
                    match = re.search(pattern, sub_line)
                    if match:
                        x = float(match.group(2))
                        y = float(match.group(3))
                        z = float(match.group(4))
                        normal_mode.append([x, y, z])
                    else:
                        break
                break

    return np.array(normal_mode), frequency


def displaced_species_along_mode(
    species: Species,
    normal_mode = np.array,
    disp_factor: float = 1.0,
    max_atom_disp: float = 99.9,
) -> Optional[Species]:
    """
    Displace the geometry along a normal mode with mode number indexed from 0,
    where 0-2 are translational normal modes, 3-5 are rotational modes and 6
    is the largest magnitude imaginary mode (if present). To displace along
    the second imaginary mode we have mode_number=7

    ---------------------------------------------------------------------------
    Arguments:
        species (autode.species.Species):
        mode_number (int): Mode number to displace along

    Keyword Arguments:
        disp_factor (float): Distance to displace (default: {1.0})

        max_atom_disp (float): Maximum displacement of any atom (Å)

    Returns:
        (autode.species.Species):

    Raises:
        (autode.exceptions.CouldNotGetProperty):
    """
    coords = species.coordinates
    disp_coords = coords.copy() + disp_factor * normal_mode

    # Ensure the maximum displacement distance any single atom is below the
    # threshold (max_atom_disp), by incrementing backwards in steps of 0.05 Å,
    # for disp_factor = 1.0 Å
    for _ in range(20):

        if (
            np.max(np.linalg.norm(coords - disp_coords, axis=1))
            < max_atom_disp
        ):
            break

        disp_coords -= (disp_factor / 20) * normal_mode

    # Create a new species from the initial
    disp_species = Species(
        name=f"{species.name}_disp",
        atoms=species.atoms.copy(),
        charge=species.charge,
        mult=species.mult,
    )
    disp_species.coordinates = disp_coords

    return disp_species


def get_ade_molecules(reactant_file, product_file, ts_guess_file, charge, multiplicity):
    """
    Load the reactant, product, and transition state molecules.

    Args:
        reactant_file (str): The name of the reactant file.
        product_file (str): The name of the product file.
        ts_guess_file (str): The name of the transition state guess file.

    Returns:
        ade.Molecule: Reactant molecule.
        ade.Molecule: Product molecule.
        ade.Molecule: Transition state molecule.
    """
    reactant = ade.Molecule(reactant_file, charge=charge, mult=multiplicity)
    product = ade.Molecule(product_file, charge=charge, mult=multiplicity)
    ts = ade.Molecule(ts_guess_file, charge=charge, mult=multiplicity)

    return reactant, product, ts


def get_xyzs(path):
    """
    Get the names of the reactant and product files.

    Args:
        path (str): The path to workdir with the xyz-files.

    Returns:
        str: The name of the reactant file.
        str: The name of the product file.
    """
    rp_path = os.path.join(path, 'rp_geometries')

    reactant_file = os.path.join(rp_path, 'reactants_geometry.xyz')
    product_file = os.path.join(rp_path, 'products_geometry.xyz') 

    return reactant_file, product_file


def get_negative_frequencies(filename, charge, solvent):
    """
    Executes an external program to calculate the negative frequencies for a given file.

    Args:
        filename (str): The name of the file to be processed.
        charge (int): The charge value for the calculation.

    Returns:
        list: A list of negative frequencies.
    """
    with open('hess.out', 'w') as out:
        if solvent is not None:
            process = subprocess.Popen(f'xtb {filename} --charge {charge} --hess --alpb {solvent}'.split(), 
                                   stderr=subprocess.DEVNULL, stdout=out)
        else:
            process = subprocess.Popen(f'xtb {filename} --charge {charge} --hess'.split(), 
                                   stderr=subprocess.DEVNULL, stdout=out)
        process.wait()
    
    neg_freq = read_negative_frequencies('g98.out')
    return neg_freq


def read_negative_frequencies(filename):
    """
    Read the negative frequencies from a file.

    Args:
        filename: The name of the file.

    Returns:
        list: The list of negative frequencies.
    """
    with open(filename, 'r') as file:
        for line in file:
            if line.strip().startswith('Frequencies --'):
                frequencies = line.strip().split()[2:]
                negative_frequencies = [freq for freq in frequencies if float(freq) < 0]
                return negative_frequencies


