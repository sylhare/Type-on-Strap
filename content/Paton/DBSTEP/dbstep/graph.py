from rdkit import Chem
from rdkit.Chem import rdmolops, rdMolDescriptors, Crippen, GraphDescriptors
import numpy as np
import pandas as pd
import pkg_resources
import sys

def crippenHContribs(mol,contribs):
    """Adds Crippen molar refractivity atomic contributions from attached H atoms to a heavy atom's contribution""" 
    res = 0.0
    ccontribs = []
    for i,at in enumerate(mol.GetAtoms()):
        if at.GetAtomicNum() != 1: #look through heavy atoms
            contrib = contribs[i]
            for nbr in at.GetNeighbors(): 
                if nbr.GetAtomicNum()==1: #look at attached hydrogens
                    contrib += contribs[nbr.GetIdx()]
            res += contrib
            ccontribs.append(contrib)
    return res,ccontribs
    
def degreeContribs(mol):
    contribs = []
    for i,at in enumerate(mol.GetAtoms()): #grab all atom contributions
        deg = at.GetDegree()
        contribs.append(deg)
    return contribs
    
def load_mcgowan():
    stream = pkg_resources.resource_stream(__name__, 'mcgowan.csv')
    return pd.read_csv(stream)

def mcgowanHContribs(molH): 
    """Adds McGowan volume atomic contributions from hydrogen to a heavy atoms contribution"""
    #Read mcgowan volumes 
    mcgowan = load_mcgowan()
    contribs = []
    for i,at in enumerate(molH.GetAtoms()):#grab all atom contributions
        num = at.GetAtomicNum()
        vol = mcgowan.loc[num]['McGowan']
        symbol = mcgowan.loc[num]['Symbol']
        contribs.append(vol)
    
    ccontribs = []
    #get heavy atom contributions
    for i,at in enumerate(molH.GetAtoms()):
        if at.GetAtomicNum() != 1: #look through heavy atoms
            contrib = contribs[i]
            for nbr in at.GetNeighbors(): 
                contrib -= 0.5 * 6.56 #subtract bonding contributions
                if nbr.GetAtomicNum()==1: #look at attached hydrogens - add to heavy atom contribs
                    contrib += contribs[nbr.GetIdx()]
                    contrib -= 0.5 * 6.56 #subtract hydrogen bond contribution
            ccontribs.append(contrib)
    return ccontribs


def make_mol_obj(line):
    toks = line.split()
    if len(toks) > 1: # expects a smiles string, followed by a property value on each line
        # parse structure from input
        smi, prop = toks[0:2]
    elif len(toks) == 1:
        smi = toks[0]
        prop=None
    mol = Chem.MolFromSmiles(smi)

    if mol is None:
        print("Warning - Could not parse SMILES for",smi,"\nEvaluating with looser constrictions.")
        mol = Chem.MolFromSmiles(smi,sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        if mol is None:
            print("Parsing Failed. Skipping this structure")
            return None,None

    return mol,prop

def mol_to_vec(input, shared_fg, voltype, max_path_length, verbose=False):
    """
	mol_to_vec function

	Allows for calculating a graph-based sum of atomic contributions in "layers"
    based off atomic connectivity ranging away from a specified common 
    functional group.

	arguments:
        input: a file listing SMILES strings or a mol object - for SMILES strings properties can optionally be provided after each string.
        shared_fg: a SMILES/SMARTS pattern shared by provided SMILES strings, or atom index to be used as a reference.
        voltype: type of volume contributions to measure, current options are "crippen", "mcgowan", or "degree".
        max_path_length: number of layers to include contributions from base functional group.
        verbose: print information to terminal.
        
    returns:
        Pandas DataFrame containing sum of atomic contributions for each layer requested in format:
            ["layer0",..."layerN","Structure" (SMILES string),"Property" (optional)]
	"""
    # computes the atomic contributions to volume/sterics at discrete number of bond lengths away from a particular atom/functional group
    mollist, y_val, vec_df, columns = [], [], [], []
    [columns.append(str(col)+'_'+voltype.lower()) for col in range(0,max_path_length)]
    if shared_fg is False:
        sys.exit("Please specify a common functional group pattern as a SMILES string using the --fg argument. Exiting.")
    
    mol_obj = []
    if isinstance(input,str):
        # create list of mol objects from smiles strings
        for line in open(input):
            mol,prop = make_mol_obj(line)
            if mol is None: continue
            mol_obj.append(mol)
            y_val.append(prop)
    else:
        # starting with mol objects, add to list
        mol_obj.append(input)
    
    for mol in mol_obj:   
        mol2vec = []  
        mat = rdmolops.GetDistanceMatrix(mol)
        # If requested indices for functional group use those instead
        try:
            base_id = int(shared_fg)
            fg_atoms = (shared_fg,)
        except ValueError:
            # the origin is defined by a particular functional group/atom of interest shared by all structures
            patt = Chem.MolFromSmarts(shared_fg)
            # if a functional group SMILES pattern is specified, the base atom is expected to be first in the smiles string
            fg_atoms = mol.GetSubstructMatch(patt)
            if len(fg_atoms) == 0:
                print("ERR: Functional group",shared_fg,"not found in molecule:",input)
                lower_fg = shared_fg.lower()
                patt = Chem.MolFromSmarts(lower_fg)
                fg_atoms = mol.GetSubstructMatch(patt)
                if len(fg_atoms) == 0:
                    print("Parsing functional group from this structure failed. Skipping this structure")
                    vec_df.append(pd.Series())
                    continue
                else:
                    if verbose:
                        print("Found by parsing functional group as",lower_fg)
        
            base_id = mol.GetSubstructMatch(patt)[0]        
        
        dist_from_base = mat[base_id]
        mollist.append(Chem.MolToSmiles(mol))

        # This uses Crippen's atomic contributions to molecular refractivity as volumes
        # Atomic contributions to logP are also available ...
        if voltype.lower() == 'crippen':
            molH = Chem.AddHs(mol)
            mrContribs = Crippen.rdMolDescriptors._CalcCrippenContribs( molH )
            logps, mrs = zip( *mrContribs )
            # condense H atom contributions to attached heavy atom
            mr, apolsCondensed = crippenHContribs(molH,mrs)
            vols = apolsCondensed
        elif voltype.lower() == 'mcgowan':
            #grab mcgowan volumes 
            molH = Chem.AddHs(mol)
            vols = mcgowanHContribs(molH)
        elif voltype.lower() == 'degree':
            vols = degreeContribs(mol)
        
        # this is the radial count up to the max_path_length
        for level in range(0,max_path_length):
            mol2vec.append(0)
            for at, dist in enumerate(dist_from_base):
                # This will try to exclude the other atoms in the defined FG apart from the base
                if dist == level and at not in fg_atoms: 
                    mol2vec[level] += vols[at]
                elif at == base_id and level == 0: # add contributions from base atom at level 0
                    """Add contributions from other base atoms here ? or not at all"""
                    mol2vec[level] += vols[at]

        # create the vector from the successive graph levels
        vec_df.append(mol2vec)
    vec_df = pd.DataFrame(vec_df,columns=columns)
    vec_df['Structure'] = np.array(mollist)
    if len(y_val) > 0:
        vec_df['Property'] = np.array(y_val)
    return vec_df
