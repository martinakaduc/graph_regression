
import os
import pickle
import pandas as pd
from tqdm import tqdm

from rdkit.Chem.rdmolfiles import SDMolSupplier

from constants import DATA_FOLDER

# Functions for reading CSVs
def readCSV(csv_file, skiprows=1):
    df = pd.read_csv(csv_file, dtype=str, skiprows=skiprows)
    return df

def findMatchingData(df, key, key_column, result_columns):
    # result_columns: list -> return dataframe
    # result_columns: string -> return series
    return df.loc[df[key_column] == key][result_columns]

def findFolder(listStr, protein, smile):
    for fol in listStr:
        if protein in fol and fol.endswith(smile):
            return fol
    return None

class ChemicalFeaturesFactory:
    """This is a singleton class for RDKit base features."""
    _instance = None

    @classmethod
    def get_instance(cls):
        try:
            from rdkit import RDConfig
            from rdkit.Chem import ChemicalFeatures
        except ModuleNotFoundError:
            raise ImportError("This class requires RDKit to be installed.")

        if not cls._instance:
            fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
            cls._instance = ChemicalFeatures.BuildFeatureFactory(fdefName)
        return cls._instance

factory = ChemicalFeaturesFactory.get_instance()

def get_feature_dict(mol):
    if mol == None:
        return {}
        
    feature_by_group = {}
    for f in factory.GetFeaturesForMol(mol):
        feature_by_group[f.GetAtomIds()] = f.GetFamily()

    feature_dict = {}
    for key in feature_by_group:
        for atom_idx in key:
            if atom_idx in feature_dict:
                feature_dict[atom_idx].append(feature_by_group[key])
            else:
                feature_dict[atom_idx] = [feature_by_group[key]]

    return feature_dict

def load_receptor_ligand_data(keys, data_agent, training=True):
    result_list = {}
    if training:
        training = "training"
        # data_agent = training_data
    else:
        training = "testing"
        # data_agent = testing_data
        
    docking_folders = os.listdir(os.path.join(DATA_FOLDER, training))
    proteins_features = {}
    
    for key in tqdm(keys):
        receptor_name, ligand_name, _ = key.split("_")
        
        # Load ligand
        ligand_smile = data_agent.ligands_smiles[ligand_name].replace('/', '-')
        docking_folder = findFolder(docking_folders, receptor_name, ligand_smile)
        assert docking_folder is not None
        ligands_sdf = SDMolSupplier("%s/%s/%s/rank1.sdf" % (DATA_FOLDER, training, docking_folder))
        ligand = ligands_sdf[0]
        ligand_feature = get_feature_dict(ligand)
        # print("Ligand %s" % ligand_name, ligand != None)
        
        # Load receptor
        # receptor = MolFromPDBFile("%s/%s.pdb" % (PDB_FOLDER, receptor_name))
        receptor = data_agent.proteins[receptor_name]
        if receptor_name not in proteins_features:
            receptor_feature = get_feature_dict(receptor)
            proteins_features[receptor_name] = receptor_feature
        else:
            receptor_feature = proteins_features[receptor_name]
        # print("Receptor %s" % receptor_name, receptor != None)
        
        result_list[key] = (ligand, receptor, ligand_feature, receptor_feature)
        
    return result_list

# Load and save
def load_and_save_data_by_keys(keys, data_agent, training=True):
    train_dict = load_receptor_ligand_data(keys, data_agent, training=training)
    if not os.path.exists("data_GMGM"):
        os.mkdir("data_GMGM")
    
    for key, data in train_dict.items():
        with open('data_GMGM/'+key, 'wb') as f:
            pickle.dump(data, f)