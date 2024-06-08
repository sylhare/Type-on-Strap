from rdkit import Chem
from itertools import combinations_with_replacement
from rdchiral.main import rdchiralRunText
from rdkit.Chem import AllChem
from math import isclose
from rdkit.Chem import rdMolTransforms
import pandas as pd
from tqdm import tqdm
import random


num_elec_dict = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Cl": 17, "Br": 35}


def unmap_unstereo_smiles(smiles):
    """Remove atom mapping and stereo-assignments from SMILES simultaneously"""
    mol = Chem.MolFromSmiles(smiles)
    [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]

    return Chem.MolToSmiles(
        Chem.MolFromSmiles(Chem.MolToSmiles(mol)), isomericSmiles=False
    )


def map_smiles(smiles):
    """Map atoms of SMILES based on index"""
    mol = Chem.MolFromSmiles(smiles)
    [atom.SetAtomMapNum(atom.GetIdx() + 1) for atom in mol.GetAtoms()]
    return Chem.MolToSmiles(mol)


def unmap_smiles(smiles):
    """Unmap atoms of SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]

    return Chem.MolToSmiles(mol)


def get_conformer(mol):
    """Get a single MMFF optimized conformer of mol-object"""
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)

        return mol.GetConformers()[0]
    except:
        return None


def get_chiral_centers_product(mol, map_num_dict, active_idx):
    """Obtain a list of map numbers, corresponding to the stereocenters of the product"""
    mol_copy = Chem.Mol(mol)
    [atom.SetAtomMapNum(0) for atom in mol_copy.GetAtoms()]
    stereocenters_map_num = []
    si = Chem.FindPotentialStereo(mol_copy)

    for element in si:
        if str(element.type) == "Atom_Tetrahedral" and element.centeredOn in active_idx:
            stereocenters_map_num.append(map_num_dict[element.centeredOn])

    return stereocenters_map_num


def get_product_isomers(stereocenters_map_num, p_mol):
    """Make combinatorial combinations of the different stereoassignments possible for the sites in the stereocenters_map_num list"""
    product_isomers = set()
    stereocenters = [
        atom.GetIdx()
        for atom in p_mol.GetAtoms()
        if atom.GetAtomMapNum() in stereocenters_map_num
    ]
    for center in stereocenters:
        p_mol.GetAtomWithIdx(center).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    product_isomers.add(Chem.MolToSmiles(p_mol))

    for l in range(len(stereocenters)):
        for subset in combinations_with_replacement(stereocenters, l):
            for center in subset:
                p_mol.GetAtomWithIdx(center).SetChiralTag(
                    Chem.ChiralType.CHI_TETRAHEDRAL_CW
                )
            product_isomers.add(Chem.MolToSmiles(p_mol))

    return product_isomers


def get_dicts_and_idx_list(p_mol, active_map_nums):
    """Get a list of active indices, i.e., atoms undergoing a change in bond order, as well as a dictionary of the map nums"""
    active_idx = [
        atom.GetIdx()
        for atom in p_mol.GetAtoms()
        if atom.GetAtomMapNum() in active_map_nums
    ]
    map_num_dict = {atom.GetIdx(): atom.GetAtomMapNum() for atom in p_mol.GetAtoms()}

    return map_num_dict, active_idx


def check_isomers(isomer_list, map_num_for_dihedrals, angle_reactants):
    """
    Generate conformers for the product isomers and compare the angle for the selected dihedral with the corresponding angle in the reactants
    In case of failure to generate isomer: assign arbitrarily high value to give lowest priority
    Select isomer if difference angle between reactants and isomer is < 15°
    If none of the isomers are sufficiently close in dihedral angle to the reactants, select the one with the lowest deviation
    """
    angle_dev_isomers = []
    for isomer in isomer_list:
        isomer_mol = Chem.MolFromSmiles(isomer)
        product_idx_for_dihedral = [
            atom.GetIdx()
            for atom in isomer_mol.GetAtoms()
            if atom.GetAtomMapNum() in map_num_for_dihedrals
        ]
        conf = get_conformer(isomer_mol)
        if conf is not None:
            angle_isomer = rdMolTransforms.GetDihedralDeg(
                conf,
                product_idx_for_dihedral[0],
                product_idx_for_dihedral[1],
                product_idx_for_dihedral[2],
                product_idx_for_dihedral[3],
            )
            # angle_reactants is either +/-0 or +/-180; approach the angle from both the positive or negative side is OK (returned angles are in range [-180,180])
            difference_angle = abs(angle_reactants) - abs(angle_isomer)
        else:
            # if angle isomer could not be determined, set the difference angle artifically high to give it lowest priority
            difference_angle = 1000
        if abs(difference_angle) < 15:
            return isomer
        else:
            angle_dev_isomers.append(abs(difference_angle))

    min_dev = min(angle_dev_isomers)
    # select the isomer with the lowest dihedral angle deviation (but make sure the isomers are sufficiently distinct in the first place)
    if min_dev < 1000 and (max(angle_dev_isomers) - min(angle_dev_isomers)) > 10:
        for idx, value in enumerate(angle_dev_isomers):
            if isclose(value, min_dev, abs_tol=1):
                return isomer_list[idx]
    else:
        raise ValueError


def get_dipolarophile_map_num_range(reactant_smiles):
    """Get the range of the atom map numbers corresponding to the dipolarophiles"""
    reactant_smiles = map_smiles(reactant_smiles)
    reactant_smiles_list = reactant_smiles.split(".")
    for smi in reactant_smiles_list:
        if "+" in smi and "-" in smi:
            continue
        else:
            map_num_list = [
                atom.GetAtomMapNum() for atom in Chem.MolFromSmiles(smi).GetAtoms()
            ]

    return range(min(map_num_list), max(map_num_list) + 1)


def get_stereo_corr_products(reactant_smiles, reaction_smarts):
    """
    Main function which generates products with RDChiral for a given reactant SMILES and SMARTS string
    and then returns a stereocompatible stereo-isomer SMILES
    """
    stereo_corr_products = []
    outcomes = rdchiralRunText(
        reaction_smarts,
        reactant_smiles,
        keep_mapnums=False,
        combine_enantiomers=True,
        return_mapped=True,
    )
    r_mol = Chem.MolFromSmiles(reactant_smiles)
    # Determine range of mapnum of dipolarophile so that you can determine whether two steroecenters are
    # in dipolarophile or not
    dipolarophile_range = get_dipolarophile_map_num_range(reactant_smiles)

    # Set atom map numbering in r_mol
    [atom.SetAtomMapNum(atom.GetIdx() + 1) for atom in r_mol.GetAtoms()]

    for product in outcomes[0]:
        # Identify all the active atoms, i.e., all those which undergo change in bonding
        active_map_nums = list(set(outcomes[1][product][1]))
        p_mol = Chem.MolFromSmiles(outcomes[1][product][0])
        map_num_dict, active_idx = get_dicts_and_idx_list(p_mol, active_map_nums)
        # Identify all the active stereocenters in the product
        stereocenters_product_map_num = get_chiral_centers_product(
            p_mol, map_num_dict, active_idx
        )

        stereocenters_from_dipolarophile_map_num = [
            map_num
            for map_num in stereocenters_product_map_num
            if map_num in dipolarophile_range
        ]
        stereocenters_dipolarophile = [
            atom.GetIdx()
            for atom in r_mol.GetAtoms()
            if atom.GetAtomMapNum() in stereocenters_from_dipolarophile_map_num
        ]

        # Get norbornene and oxo-norbornadiene patterns since stereochemistry for these types of reactions does not need to be constrained
        mol_norbornene = Chem.MolFromSmiles("C1=CC2CCC1C2")
        mol_oxo_norbornadiene = Chem.MolFromSmiles("C1=CC2C=CC1O2")

        # If only a single stereocenter is present in the dipolarophile, or if dipolarophile is rigid (norbornene and oxo-norbornadiene),
        # then this doesn't matter and no further checks are needed, lowest product conformer can be taken without restrictions
        if (
            len(stereocenters_from_dipolarophile_map_num) <= 1
            or len(r_mol.GetSubstructMatch(mol_norbornene)) != 0
            or len(r_mol.GetSubstructMatch(mol_oxo_norbornadiene)) != 0
        ):
            stereo_corr_products.append(Chem.MolToSmiles(p_mol))
        # If two stereocenters in dipolarophile, select compatible isomer (all other stereocenters can be set independently by autodE)
        elif len(stereocenters_from_dipolarophile_map_num) == 2:
            conf_r = get_conformer(r_mol)
            neighbors_cistrans = []
            # First try to select a neighbor that is not in ring (both neighbors in ring is problematic for dihedral determination),
            # if this doesn't work, then settle for a ring-atom
            for center in stereocenters_dipolarophile:
                try:
                    neighbors_cistrans.append(
                        [
                            atom
                            for atom in r_mol.GetAtomWithIdx(center).GetNeighbors()
                            if (
                                atom.GetIdx() not in stereocenters_dipolarophile
                                and atom.IsInRing() == False
                            )
                        ][0]
                    )
                except IndexError:
                    neighbors_cistrans.append(
                        [
                            atom
                            for atom in r_mol.GetAtomWithIdx(center).GetNeighbors()
                            if atom.GetIdx() not in stereocenters_dipolarophile
                        ][0]
                    )

            map_num_for_dihedral = list(
                map(
                    lambda x: x.GetAtomMapNum(),
                    [
                        neighbors_cistrans[0],
                        r_mol.GetAtomWithIdx(stereocenters_dipolarophile[0]),
                        r_mol.GetAtomWithIdx(stereocenters_dipolarophile[1]),
                        neighbors_cistrans[1],
                    ],
                )
            )

            reactant_idx_for_dihedral = [
                atom.GetIdx()
                for atom in r_mol.GetAtoms()
                if atom.GetAtomMapNum() in map_num_for_dihedral
            ]
            dihedral_r = rdMolTransforms.GetDihedralDeg(
                conf_r,
                reactant_idx_for_dihedral[0],
                reactant_idx_for_dihedral[1],
                reactant_idx_for_dihedral[2],
                reactant_idx_for_dihedral[3],
            )

            product_isomers = get_product_isomers(
                stereocenters_from_dipolarophile_map_num, p_mol
            )

            try:
                stereo_corr_product = check_isomers(
                    list(product_isomers), map_num_for_dihedral, dihedral_r
                )
                stereo_corr_products.append(stereo_corr_product)
            except ValueError:
                print(f"Error for {reactant_smiles}")
                continue

    return stereo_corr_products


def get_rxn_smiles(reactant_smiles, smarts_list, rxn_smiles_list):
    """
    Iterate through the SMARTS list for a given reactant SMILES,
    try to get stereo compatible product and make sure the resulting (unmapped, unstereo'ed) SMILES
    has not previously been generated
    """
    unmapped_products_no_stereo = []
    for smarts in smarts_list:
        try:
            corr_prods = get_stereo_corr_products(reactant_smiles, smarts)
            for corr_prod in corr_prods:
                if unmap_unstereo_smiles(corr_prod) not in unmapped_products_no_stereo:
                    rxn_smiles_list.append(
                        f"{map_smiles(reactant_smiles)}>>{corr_prod}"
                    )
                    unmapped_products_no_stereo.append(unmap_unstereo_smiles(corr_prod))
                else:
                    continue
        except IndexError:
            continue

    return rxn_smiles_list


def get_data_point(idx, rxn_smiles, solvent, temp):
    """Return data point info in correct format for writing to file"""
    return f"{idx},{rxn_smiles},{solvent},{temp} \n"


def get_bio_dipolarophiles_list(n=1500):
    """Sample the list of biofragment-based dipolarophiles"""
    bio_dipolarophiles_list = [
        "C=C",
        "C/C=C/C",
        "C/C=C\C",  # prostaglandines etc.
        "[O-]C(=O)/C=C\C(=O)[O-]",
        "[O-]C(=O)/C=C/C(=O)[O-]",  # fumaric and maleic acid
        "C1(=O)C=CC(=O)C=C1",
        "C1(=O)C(OC)=C(OC)C(=O)C(C)=C1C",  # ubiquinone
        "CC(C)=C(C)C",  # terpineol
        "C1C=CN(C)C=C1C(N)=O",  # NADH
        "CC(C)=CC=O",  # retinal
        "CC(C)=C/C=C/C(C)=C",  # retinol, beta-caroteen
        "[NH2+]=C(N)NC",  # argininosuccinic acid
    ]

    df_biofragments = pd.DataFrame(bio_dipolarophiles_list)
    df_biofragments_sample = df_biofragments.sample(n=n, replace=True, random_state=3)
    biofragments_list = df_biofragments_sample[0].values.tolist()

    return biofragments_list


def get_dipole_list(filename="sample_lists/dipoles_sample.csv"):
    """Get the dipole list from input file"""
    dipoles = pd.read_csv(filename)
    dipoles_list = dipoles["dipole"].values.tolist()

    return dipoles_list


def get_synthetic_dipolarophiles_list(
    filename="sample_lists/dipolarophiles_sample.csv",
):
    """Get the (shuffled) dipolarophile list from input file"""
    dipolarophiles = pd.read_csv(filename)
    dipolarophiles = dipolarophiles.sample(frac=1, random_state=3)
    dipolarophiles_list = dipolarophiles["dipolarophile"].values.tolist()

    return dipolarophiles_list


def get_smarts_lists():
    """Get the SMARTS strings for each of the main types of dipole/dipolarophiles patterns"""
    (
        smarts_list_double,
        smarts_list_triple,
        smarts_list_triple2,
        smarts_list_aromatic,
    ) = ([], [], [], [])
    smarts_list_triple.append(
        "([*:1]#[*+:2][*-:3].[C,N:4]=[C:5])>>[*:1]1=[*;+0:2][*;-0:3][C,N:4][C:5]1"
    )
    smarts_list_triple.append(
        "([*:1]#[*+:2][*-:3].[C,N:4]=[C:5])>>[*:1]1=[*;+0:2][*;-0:3][C:5][C,N:4]1"
    )
    smarts_list_triple.append(
        "([*:1]=[*+:2][*-:3].[C:4]#[C:5])>>[*:1]1[*;+0:2][*;-0:3][C:4]=[C:5]1"
    )
    smarts_list_triple.append(
        "([*:1]=[*+:2][*-:3].[C:4]#[C:5])>>[*:1]1[*;+0:2][*;-0:3][C:5]=[C:4]1"
    )

    smarts_list_triple2.append(
        "([*:1]#[*+:2][*-:3].[C:4]#[C:5])>>[*:1]1=[*;+0:2][*;-0:3][C:4]=[C:5]1"
    )
    smarts_list_triple2.append(
        "([*:1]#[*+:2][*-:3].[C:4]#[C:5])>>[*:1]1=[*;+0:2][*;-0:3][C:5]=[C:4]1"
    )

    smarts_list_double.append(
        "([*:1]=[*+:2][*-:3].[C,N:4]=[C:5])>>[*:1]1[*;+0:2][*;-0:3][C,N:4][C:5]1"
    )
    smarts_list_double.append(
        "([*:1]=[*+:2][*-:3].[C,N:4]=[C:5])>>[*:1]1[*;+0:2][*;-0:3][C:5][C,N:4]1"
    )

    smarts_list_aromatic.append(
        "([*:1]:[*+:2][*-:3].[C,N:4]=[C:5])>>[*:1]1[*;+0:2][*;-0:3][C,N:4][C:5]1"
    )
    smarts_list_aromatic.append(
        "([*:1]:[*+:2][*-:3].[C,N:4]=[C:5])>>[*:1]1[*;+0:2][*;-0:3][C:5][C,N:4]1"
    )
    smarts_list_aromatic.append(
        "([*:1]:[*+:2][*-:3].[C:4]#[C:5])>>[*:1]1[*;+0:2][*;-0:3][C:4]=[C:5]1"
    )
    smarts_list_aromatic.append(
        "([*:1]:[*+:2][*-:3].[C:4]#[C:5])>>[*:1]1[*;+0:2][*;-0:3][C:5]=[C:4]1"
    )

    return (
        smarts_list_double,
        smarts_list_triple,
        smarts_list_triple2,
        smarts_list_aromatic,
    )


def get_r_smiles_list(dipole_list, dipolarophiles_list):
    """Combine dipole and dipolarophile lists into reactant SMILES list"""
    r_smiles_list = []
    for idx, dipole in enumerate(dipole_list):
        r_smiles_list.append(f"{dipolarophiles_list[idx]}.{dipole}")

    return r_smiles_list


def get_num_elec(reactant_smiles):
    """Get the number of electrons in the reactant system based on dictionary"""
    mol = Chem.MolFromSmiles(reactant_smiles)
    mol = Chem.AddHs(mol)

    return sum([num_elec_dict[atom.GetSymbol()] for atom in mol.GetAtoms()])


def reorder_entries_based_on_num_elec(r_smiles_list):
    """Reorder the entries based on their size to make sure each batch of reaction profile computation takes a similar amount of time"""
    df = pd.DataFrame(r_smiles_list)
    df["reactant_smiles"] = df[0]
    df["num_elec"] = df["reactant_smiles"].apply(
        lambda x: get_num_elec(x.split(">")[0])
    )
    df.sort_values("num_elec", inplace=True)
    r_smiles_list = df["reactant_smiles"].tolist()

    return r_smiles_list


def reorder_reactants(reactant_smiles):
    """Reorder the reactants in a reactant SMILES so that the dipolarophile always comes first"""
    reactant_smiles_list = reactant_smiles.split(".")
    if "+" not in reactant_smiles_list[0] and "-" not in reactant_smiles_list[0]:
        return f"{reactant_smiles_list[0]}.{reactant_smiles_list[1]}"
    else:
        return f"{reactant_smiles_list[1]}.{reactant_smiles_list[0]}"


if __name__ == "__main__":
    # initialization for synthetic reactions
    rxn_smiles_list = []
    dipole_list = get_dipole_list("work_dir/dipoles.csv")
    synthetic_dipolarophiles_list = get_synthetic_dipolarophiles_list(
        "work_dir/dipolarophiles.csv"
    )
    (
        smarts_list_double,
        smarts_list_triple,
        smarts_list_triple2,
        smarts_list_aromatic,
    ) = get_smarts_lists()

    r_smiles_list = get_r_smiles_list(random.sample(dipole_list, 30), random.sample(synthetic_dipolarophiles_list, 30))
    r_smiles_list = reorder_entries_based_on_num_elec(list(set(r_smiles_list)))

    # get the reaction SMILES for each reactant combination
    for r_smiles in tqdm(r_smiles_list):
        if "n" in r_smiles or "o" in r_smiles:
            rxn_smiles_list = get_rxn_smiles(
                r_smiles, smarts_list_aromatic, rxn_smiles_list
            )
        elif r_smiles.count("#") - r_smiles.count("C#N") - r_smiles.count("N#C") == 0:
            rxn_smiles_list = get_rxn_smiles(
                r_smiles, smarts_list_double, rxn_smiles_list
            )
        elif r_smiles.count("#") - r_smiles.count("C#N") - r_smiles.count("N#C") == 1:
            rxn_smiles_list = get_rxn_smiles(
                r_smiles, smarts_list_triple, rxn_smiles_list
            )
        elif r_smiles.count("#") - r_smiles.count("C#N") - r_smiles.count("N#C") == 2:
            rxn_smiles_list = get_rxn_smiles(
                r_smiles, smarts_list_triple2, rxn_smiles_list
            )

    # write to file
    rxn_smiles_written_synthetic, i = [], 0
    with open("work_dir/reactions_synthetic_finalized.csv", "w") as f:
        f.write(",rxn_smiles,solvent,temp\n")
        for rxn_smiles in rxn_smiles_list:
            if rxn_smiles not in rxn_smiles_written_synthetic:
                f.write(get_data_point(i, rxn_smiles, "water", 298.15))
                rxn_smiles_written_synthetic.append(rxn_smiles)
                i += 1
            else:
                continue
