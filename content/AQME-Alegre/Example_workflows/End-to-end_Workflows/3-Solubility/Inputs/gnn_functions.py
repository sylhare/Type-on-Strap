import tensorflow as tf
from tensorflow.keras import layers
import nfp
from nfp.preprocessing.xtb_preprocessor import xTBSmilesPreprocessor
from nfp.preprocessing.features import get_ring_size
import numpy as np
from tqdm import tqdm
tqdm.pandas()


def atom_featurizer(atom):
    """ Return a string representing the atom type
    """

    return str((
        atom.GetSymbol(),
        atom.GetIsAromatic(),
        get_ring_size(atom, max_size=6),
        atom.GetDegree(),
        atom.GetTotalNumHs(includeNeighbors=True)
    ))

def bond_featurizer(bond, flipped=False):
    """ Get a similar classification of the bond type.
    Flipped indicates which 'direction' the bond edge is pointing. """
    
    if not flipped:
        atoms = "{}-{}".format(
            *tuple((bond.GetBeginAtom().GetSymbol(),
                    bond.GetEndAtom().GetSymbol())))
    else:
        atoms = "{}-{}".format(
            *tuple((bond.GetEndAtom().GetSymbol(),
                    bond.GetBeginAtom().GetSymbol())))
    
    btype = str(bond.GetBondType())
    ring = 'R{}'.format(get_ring_size(bond, max_size=6)) if bond.IsInRing() else ''
    
    return " ".join([atoms, btype, ring]).strip()

def data_split(valid, test, train, preprocessor, sol):
    train_dataset = tf.data.Dataset.from_generator(
        lambda: ((preprocessor(row.smiles, row.xtbjson, train=True), row['measured log solubility in mols per litre'])
                for i, row in sol[sol.smiles.isin(train.smiles)].iterrows()),
        output_signature=(preprocessor.output_signature, tf.TensorSpec((), dtype=tf.float32)))\
        .cache().shuffle(buffer_size=200)\
        .padded_batch(batch_size=(len(train)), padding_values=(preprocessor.padding_values, tf.constant(np.nan, dtype=tf.float32)))\
        .prefetch(tf.data.experimental.AUTOTUNE)

    valid_dataset = tf.data.Dataset.from_generator(
        lambda: ((preprocessor(row.smiles, row.xtbjson, train=False), row['measured log solubility in mols per litre'])
                for i, row in sol[sol.smiles.isin(valid.smiles)].iterrows()),
        output_signature=(preprocessor.output_signature, tf.TensorSpec((), dtype=tf.float32)))\
        .cache()\
        .padded_batch(batch_size=(len(valid)))\
        .prefetch(tf.data.experimental.AUTOTUNE)

    test_dataset = tf.data.Dataset.from_generator(
        lambda: (preprocessor(row.smiles, row.xtbjson, train=False)
                for i, row in test.iterrows()),
        output_signature=preprocessor.output_signature)\
        .padded_batch(batch_size=len(test))\
        .prefetch(tf.data.experimental.AUTOTUNE)
    
    return train_dataset, valid_dataset, test_dataset

xtb_mol_features = ['HOMO','LUMO','Dipole module/D','Total charge','Fermi-level/eV','Total dispersion C6','Total dispersion C8','Total polarizability alpha']
# Define how to featurize the input molecules
preprocessor = xTBSmilesPreprocessor(atom_features=atom_featurizer, bond_features=bond_featurizer, xtb_bond_features=[], xtb_mol_features=xtb_mol_features, scaler=False)

def gnn_model():

    # Define the keras model
    atom = layers.Input(shape=[None], dtype=tf.int64, name="atom")
    bond = layers.Input(shape=[None], dtype=tf.int64, name="bond")
    connectivity = layers.Input(shape=[None, 2], dtype=tf.int64, name="connectivity")
    atom_xtb = layers.Input(shape=[None, len(preprocessor.xtb_atom_features)], dtype=tf.float64, name="atom_xtb")
    mol_xtb = layers.Input(shape=[len(preprocessor.xtb_mol_features), ], dtype=tf.float64, name="mol_xtb")
    atom_features = 128

    # Initialize the atom states
    atom_class = layers.Embedding(
        preprocessor.atom_classes,
        atom_features, name='atom_embedding', mask_zero=True)(atom)

    atom_class_xtb = layers.Dense(atom_features,activation='relu')(atom_xtb)
    atom_state = layers.Concatenate(axis=-1)([atom_class, atom_class_xtb])
    atom_state = layers.Dense(atom_features)(atom_state)

    # Initialize the bond states
    bond_state = layers.Embedding(
        preprocessor.bond_classes,
        atom_features, name='bond_embedding', mask_zero=True)(bond)
    
    mol_state = layers.Dense(atom_features,activation='relu')(mol_xtb)

    # Here we use our first nfp layer. This is an attention layer that looks at
    # the atom and bond states and reduces them to a single, graph-level vector. 
    # mum_heads * units has to be the same dimension as the atom / bond dimension
    global_state = nfp.GlobalUpdate(units=atom_features, num_heads=1)([atom_state, bond_state, connectivity, mol_state])

    for _ in range(6):  # Do the message passing
        new_bond_state = nfp.EdgeUpdate()([atom_state, bond_state, connectivity, global_state])
        bond_state = layers.Add()([bond_state, new_bond_state])
        
        new_atom_state = nfp.NodeUpdate()([atom_state, bond_state, connectivity, global_state])
        atom_state = layers.Add()([atom_state, new_atom_state])
        
        new_global_state = nfp.GlobalUpdate(units=atom_features, num_heads=1)(
            [atom_state, bond_state, connectivity, global_state]) 
        global_state = layers.Add()([global_state, new_global_state])
    
    # Since the final prediction is a single, molecule-level property (solubility), we 
    # reduce the last global state to a single prediction.
    sol_prediction = layers.Dense(1)(global_state)

    # Construct the tf.keras model
    model = tf.keras.Model([atom, bond, connectivity, atom_xtb, mol_xtb], [sol_prediction])
    
    return model

def gnn_data(valid, test, train, sol):
    
    # Define train, validation and test datasets
    train_dataset, valid_dataset, test_dataset = data_split(valid, test, train, preprocessor, sol)
    
    return train_dataset, valid_dataset, test_dataset