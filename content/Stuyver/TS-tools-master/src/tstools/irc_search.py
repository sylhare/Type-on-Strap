import re
import autode as ade
import os
from autode.mol_graphs import make_graph
import subprocess

atomic_number_to_symbol = {
    1: 'H',  2: 'He',  3: 'Li',  4: 'Be',  5: 'B',
    6: 'C',  7: 'N',   8: 'O',   9: 'F',  10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
    16: 'S',  17: 'Cl', 18: 'Ar', 19: 'K',  20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V',  24: 'Cr', 25: 'Mn',
    26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br',
    36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y',  40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh',
    46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I',  54: 'Xe', 55: 'Cs',
    56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb',
    66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W',  75: 'Re',
    76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
    86: 'Rn'
}


def extract_irc_geometries(log_path_forward, log_path_reverse):
    """
    Extract IRC (Intrinsic Reaction Coordinate) geometries from forward and reverse log files,
    and write them to XYZ files.

    Parameters:
    - log_path_forward (str): Path to the forward IRC log file.
    - log_path_reverse (str): Path to the reverse IRC log file.
    """
    geometry_block_forward = extract_geometry_block_from_irc(log_path_forward)
    geometry_block_reverse = extract_geometry_block_from_irc(log_path_reverse)

    write_geometry_block_to_xyz(geometry_block_forward, f'{log_path_forward[:-4]}.xyz', True)
    write_geometry_block_to_xyz(geometry_block_reverse, f'{log_path_reverse[:-4]}.xyz', True)


def extract_geometry_block_from_irc(log_path):
    """
    Extract the geometry block from an IRC (Intrinsic Reaction Coordinate) log file.

    Parameters:
    - log_path (str): Path to the IRC log file.

    Returns:
    List[str]: List of lines containing the geometry block.
    """
    geometry_start_pattern = re.compile(r'^\s*Cartesian Coordinates \(Ang\):\s*$')
    geometry_end_pattern = re.compile(r'^\s*CHANGE IN THE REACTION COORDINATE =\s*([+-]?\d*\.\d+)\s*$')
    with open(log_path, 'r') as log_file:
        lines = log_file.readlines()

        for i, line in enumerate(lines):
            if geometry_start_pattern.match(line.strip()):
                geometry_block_start_line = i + 5
            if geometry_end_pattern.search(line):
                geometry_block_end_line = i - 1

    return lines[geometry_block_start_line:geometry_block_end_line] 


def write_geometry_block_to_xyz(geometry_block, output_xyz_path, irc=False):
    """
    Write a geometry block to an XYZ file.

    Parameters:
    - geometry_block (List[str]): List of lines containing the geometry block.
    - output_xyz_path (str): Path to the output XYZ file.
    - irc (bool): Indicates whether the geometry block is from an IRC log. Defaults to False.
    """
    with open(output_xyz_path, 'w') as xyz_file:
        # Write the number of atoms
        xyz_file.write(str(len(geometry_block)) + '\n\n')
        # Write the atomic coordinates
        for line in geometry_block:
            split_line = line.split()
            if irc:
                xyz_file.write(f'{atomic_number_to_symbol[int(split_line[1])]} {float(split_line[2]):.6f} {float(split_line[3]):.6f} {float(split_line[4]):.6f}\n')
            else:
                xyz_file.write(f'{atomic_number_to_symbol[int(split_line[1])]} {float(split_line[3]):.6f} {float(split_line[4]):.6f} {float(split_line[5]):.6f}\n')
    

def extract_transition_state_geometry(log_path, output_xyz_path):
    """
    Extract the geometry of a transition state from a Gaussian16 log file and save it to an XYZ file.

    Parameters:
    - log_path (str): Path to the Gaussian16 log file.
    - output_xyz_path (str): Path to the output XYZ file.
    """
    # Define regular expressions to identify relevant lines
    geometry_start_pattern = re.compile(r'^\s*Standard orientation:.*')
    geometry_end_pattern = re.compile(r'^\s*-+\s*$')
    transition_state_pattern = re.compile(r'^\s*-- Stationary point found\.$')

    # Set boolean flag
    ts_found = False
    geometry_block_start_line = 0
    geometry_block_end_line = 0

    # Read the Gaussian16 log file
    with open(log_path, 'r') as log_file:
        lines = log_file.readlines()

        # Iterate through the lines in the log file
        for i, line in enumerate(lines):
            # Break if transition state is found
            if transition_state_pattern.match(line):
                ts_found = True

            if ts_found and geometry_start_pattern.match(line):
                geometry_block_start_line = i + 5
            
            if ts_found and geometry_block_start_line != 0 and i > geometry_block_start_line and geometry_end_pattern.match(line):
                geometry_block_end_line = i
                break

    geometry_block = lines[geometry_block_start_line: geometry_block_end_line]

    # Write the final geometry to an XYZ file
    write_geometry_block_to_xyz(geometry_block, f'{log_path[:-4]}.xyz')

    print(f'Transition state geometry extracted and saved to {output_xyz_path}.')


def write_coordinates_to_xyz(coordinates, file_path):
    """
    Write atomic coordinates to an XYZ file.

    Parameters:
    - coordinates (list): A list of tuples representing atomic coordinates.
        Each tuple should have the format (atom_symbol, x, y, z).
    - file_path (str): The path to the output XYZ file.

    Returns:
    None
    """
    with open(file_path, 'w') as xyz_file:
        # Write the number of atoms
        xyz_file.write(str(len(coordinates)) + '\n\n')
        for coordinate in coordinates:
            xyz_file.write(f'{coordinate[0]} {float(coordinate[1]):.6f} {float(coordinate[2]):.6f} {float(coordinate[3]):.6f}\n')


def extract_coordinates(file_path):
    """
    Extract atomic coordinates from a file and write them to a new XYZ file.

    Parameters:
    - file_path (str): Path to the input file.

    Returns:
    None
    """
    coordinates = []

    with open(file_path, 'r') as file:
        found_coordinates = False

        for line in file:
            if 'final structure' in line.lower():
                found_coordinates = True
                continue

            if found_coordinates and line.strip() == '':
                break

            if found_coordinates:
                atom_info = line.split()
                if len(atom_info) == 4:
                    element, x, y, z = atom_info
                    coordinates.append((element, float(x), float(y), float(z)))

    write_coordinates_to_xyz(coordinates, f'{file_path[:-4]}_opt.xyz')


def write_xtb_input_file(xyz_path):
    """
    Write an XTB input file based on an XYZ file.

    Parameters:
    - xyz_path (str): Path to the XYZ file.

    Returns:
    str: Path to the created XTB input file.
    """
    file_path = os.path.join(os.path.dirname(xyz_path), 'xtb.inp')
    input_block = """
    $wall
    potential=logfermi
    sphere: auto, all
    $end
    """

    with open(file_path, 'w') as file:
        file.write(input_block)

    return file_path


def optimize_final_point_irc(xyz_file, charge, multiplicity, solvent=None):
    """
    Optimize the final point of an Intrinsic Reaction Coordinate (IRC) using XTB.

    Parameters:
    - xyz_file (str): Path to the XYZ file containing the molecular coordinates.
    - charge (int): Charge of the system.
    - multiplicity (int): Spin multiplicity of the system (1 for singlet, 2 for doublet, etc.).
    - solvent (str, optional): Solvent information. Defaults to None.

    Returns:
    None: The function writes the optimized coordinates to an output file.
    """
    inp_path = write_xtb_input_file(xyz_file)
    with open(f'{xyz_file[:-4]}.out', 'w') as out:
        cmd = f'xtb --input {inp_path} {xyz_file} --opt --cma --charge {charge} '

        if multiplicity == 1:
            pass
        elif multiplicity == 2:
            cmd += '--uhf 1 '
        
        if solvent is not None:
            cmd += f'--alpb {solvent} '

        process = subprocess.Popen(cmd.split(), stderr=subprocess.DEVNULL, stdout=out)
        process.wait()

    extract_coordinates(f'{xyz_file[:-4]}.out')


def update_molecular_graphs(rel_tolerance, forward_mol, reverse_mol, reactant_mol, product_mol):
    """
    Update molecular graphs for a set of molecules using a relative tolerance.

    Parameters:
    - rel_tolerance (float): Relative tolerance for updating molecular graphs.
    - forward_mol (ade.Molecule): Forward molecule.
    - reverse_mol (ade.Molecule): Reverse molecule.
    - reactant_mol (ade.Molecule): Reactant molecule.
    - product_mol (ade.Molecule): Product molecule.

    Returns:
    Tuple[ade.Molecule, ade.Molecule, ade.Molecule, ade.Molecule]: Updated forward, reverse, reactant, and product molecules.
    """
    for mol in [forward_mol, reverse_mol, reactant_mol, product_mol]:
        make_graph(mol, rel_tolerance=rel_tolerance)

    return forward_mol, reverse_mol, reactant_mol, product_mol 


#TODO: optimize in case of DFT???
def compare_molecules_irc(forward_xyz, reverse_xyz, reactant_xyz, product_xyz, charge=0, multiplicity=1, solvent=None):
    """
    Compare molecular graphs of IRC (Intrinsic Reaction Coordinate) molecules after optimization.

    Parameters:
    - forward_xyz (str): Path to the XYZ file for the forward IRC molecule.
    - reverse_xyz (str): Path to the XYZ file for the reverse IRC molecule.
    - reactant_xyz (str): Path to the XYZ file for the reactant molecule.
    - product_xyz (str): Path to the XYZ file for the product molecule.
    - charge (int, optional): Charge of the system. Defaults to 0.
    - multiplicity (int, optional): Spin multiplicity of the system (1 for singlet, 2 for doublet, etc.). Defaults to 1.
    - solvent (str, optional): Solvent information. Defaults to None.

    Returns:
    bool: True if the molecular graphs match for any relative tolerance, False otherwise.
    """
    # first reoptimize the final points
    optimize_final_point_irc(forward_xyz, charge, multiplicity, solvent)
    optimize_final_point_irc(reverse_xyz, charge, multiplicity, solvent)
    # then take final geometry and do actual comparison
    forward_mol = ade.Molecule(f'{forward_xyz[:-4]}_opt.xyz', name='forward', charge=charge, mult=multiplicity)
    reverse_mol = ade.Molecule(f'{reverse_xyz[:-4]}_opt.xyz', name='reverse', charge=charge, mult=multiplicity)
    reactant_mol = ade.Molecule(reactant_xyz, name='reactant', charge=charge, mult=multiplicity)
    product_mol = ade.Molecule(product_xyz, name='product', charge=charge, mult=multiplicity)

    # if no bond change between reactant and product, abort
    if set(reactant_mol.graph.edges) == set(product_mol.graph.edges):
        return False

    for rel_tolerance in [0.3, 0.25, 0.20, 0.15, 0.10, 0.05, 0.01]:
        forward_mol, reverse_mol, reactant_mol, product_mol = update_molecular_graphs(rel_tolerance, forward_mol, reverse_mol, reactant_mol, product_mol)
        if set(forward_mol.graph.edges) == set(reactant_mol.graph.edges) and set(reverse_mol.graph.edges) == set(product_mol.graph.edges) \
            and set(forward_mol.graph.edges) != set(reverse_mol.graph.edges):
            return True
        elif set(forward_mol.graph.edges) == set(product_mol.graph.edges) and set(reverse_mol.graph.edges) == set(reactant_mol.graph.edges) \
            and set(forward_mol.graph.edges) != set(reverse_mol.graph.edges): 
            return True
        else:
            continue
    
    return False


def generate_gaussian_irc_input(xyz_file, output_prefix='irc_calc', method='B3LYP/6-31G**', 
                                mem='2GB', proc=2, solvent=None, charge=0, multiplicity=1):
    """
    Generate Gaussian input files for IRC (Intrinsic Reaction Coordinate) calculations.

    Parameters:
    - xyz_file (str): Path to the XYZ file containing atomic coordinates.
    - output_prefix (str, optional): Prefix for output files. Defaults to 'irc_calc'.
    - method (str, optional): Gaussian method and basis set. Defaults to 'B3LYP/6-31G**'.
    - mem (str, optional): Gaussian memory specification. Defaults to '2GB'.
    - proc (int, optional): Number of processors. Defaults to 2.
    - solvent (str, optional): Solvent information. Defaults to None.
    - charge (int, optional): Charge of the system. Defaults to 0.
    - multiplicity (int, optional): Spin multiplicity of the system. Defaults to 1.

    Returns:
    Tuple[str, str]: Paths to the forward and reverse Gaussian input files.
    """
    # Read the XYZ file
    with open(xyz_file, 'r') as xyz:
        lines = xyz.readlines()

    # Extract the atomic coordinates
    atom_coords = lines[2:]

    # Define the Gaussian input file content
    if 'external' not in method:
        if solvent is not None:
            input_content_f = f'%Chk={xyz_file.split("/")[-1][:-4]}.chk\n%NProc={proc}\n%Mem={mem}\n#p IRC(calcfc, maxpoint=50, stepsize=15, Forward) {method} SCRF=(Solvent={solvent}, smd)' \
                f'\n\nIRC Calculation\n\n{charge} {multiplicity}\n{"".join(atom_coords)}\n\n'
            input_content_r = f'%Chk={xyz_file.split("/")[-1][:-4]}.chk\n%NProc={proc}\n%Mem={mem}\n#p IRC(calcfc, maxpoint=50, stepsize=15, Reverse) {method} SCRF=(Solvent={solvent}, smd)' \
                f'\n\nIRC Calculation\n\n{charge} {multiplicity}\n{"".join(atom_coords)}\n\n' 
        else:
            input_content_f = f'%Chk={xyz_file.split("/")[-1][:-4]}.chk\n%NProc={proc}\n%Mem={mem}\n#p IRC(calcfc, maxpoint=50, stepsize=15, Forward) {method}' \
                f'\n\nIRC Calculation\n\n{charge} {multiplicity}\n{"".join(atom_coords)}\n\n'
            input_content_r = f'%Chk={xyz_file.split("/")[-1][:-4]}.chk\n%NProc={proc}\n%Mem={mem}\n#p IRC(calcfc, maxpoint=50, stepsize=15, Reverse) {method}' \
                f'\n\nIRC Calculation\n\n{charge} {multiplicity}\n{"".join(atom_coords)}\n\n' 
    else:
        if solvent is not None:
            external_method = f'external="{method} alpb={solvent}"'
        else:
            external_method = f'external="{method}"'
        input_content_f = f'%Chk={xyz_file.split("/")[-1][:-4]}.chk\n#p IRC(calcfc, maxpoint=50, stepsize=15, Forward) {external_method}\n\n' \
            f'IRC Calculation\n\n{charge} {multiplicity}\n{"".join(atom_coords)}\n\n'
        input_content_r = f'%Chk={xyz_file.split("/")[-1][:-4]}.chk\n#p IRC(calcfc, maxpoint=50, stepsize=15, Reverse) {external_method}\n\n' \
            f'IRC Calculation\n\n{charge} {multiplicity}\n{"".join(atom_coords)}\n\n' 

    # Write the input content to a Gaussian input file -- forward
    input_filename_f = os.path.join(os.path.dirname(xyz_file), f'{output_prefix}_forward.com')
    with open(input_filename_f, 'w') as input_file:
        input_file.write(input_content_f)
    
    # Write the input content to a Gaussian input file -- reverse
    input_filename_r = os.path.join(os.path.dirname(xyz_file), f'{output_prefix}_reverse.com')
    with open(input_filename_r, 'w') as input_file:
        input_file.write(input_content_r)

    return input_filename_f, input_filename_r

if __name__ == '__main__':
    #path = '/Users/thijsstuyver/Desktop/reaction_profile_generator/benchmarking_2.0_200/final_ts_guesses/reaction_R100_final_ts_guess_0.xyz'
    #path = '/Users/thijsstuyver/Desktop/reaction_profile_generator/lol.log'
    #extract_transition_state_geometry(path, f'{path[:-4]}.xyz')
    #generate_gaussian_irc_input(f'{path[:-4]}.xyz', method='external="/home/thijs/Jensen_xtb_gaussian/profiles_test/extra/xtb_external.py"')
    #extract_irc_geometries('/Users/thijsstuyver/Desktop/reaction_profile_generator/lol/test_irc_forward.log', 
    #                       '/Users/thijsstuyver/Desktop/reaction_profile_generator/lol/test_irc_forward.log')
    # optimize_final_point_irc('lol/ts_guess_4_irc_forward.xyz', 0)
    reaction_correct = compare_molecules_irc(
        '/Users/thijsstuyver/Desktop/reaction_R8/g16_dir/ts_guess_2_irc_forward.xyz', 
        '/Users/thijsstuyver/Desktop/reaction_R8/g16_dir/ts_guess_2_irc_reverse.xyz',
        '/Users/thijsstuyver/Desktop/reaction_R8/rp_geometries/reactants_geometry.xyz', 
        '/Users/thijsstuyver/Desktop/reaction_R8/rp_geometries/products_geometry.xyz')
    print(reaction_correct)
    #print(reaction_correct)
    #extract_transition_state_geometry('logs/ts_guess_0.log', 'logs/test.xyz')
