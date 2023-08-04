import os
import sys
import pickle

from constants import DATA_FOLDER
from datasets import IUPHAR
from protein_utils import getProteinsByFamilies, getTrainTestByProteins
from utils import load_and_save_data_by_keys

if __name__ == '__main__':
    after_docking = sys.argv[1] == 'True'
    load_saved_pickles = sys.argv[2] == 'True'
    
    if load_saved_pickles:
        print('Loading saved pickles...')
        
        with open(os.path.join(DATA_FOLDER, 'training_proteins.pkl'), 'rb') as f:
            TRAINING_PROTEINS = pickle.load(f)
            
        with open(os.path.join(DATA_FOLDER, 'testing_proteins.pkl'), 'rb') as f:
            TESTING_PROTEINS = pickle.load(f)
            
        with open(os.path.join(DATA_FOLDER, 'proteins_by_families.pkl'), 'rb') as f:
            proteins_by_families = pickle.load(f)
            
        with open(os.path.join(DATA_FOLDER, 'protein_family_mapping.pkl'), 'rb') as f:
            protein_family_mapping = pickle.load(f)
            
        with open(os.path.join(DATA_FOLDER, 'uniprot2pdb_mapping_dict.pkl'), 'rb') as f:
            uniprot2pdb_mapping_dict = pickle.load(f)
            
    else:    
        proteins_by_families, protein_family_mapping, uniprot2pdb_mapping_dict = \
            getProteinsByFamilies(os.path.join(DATA_FOLDER, "targets_and_families.csv"))
            
        TRAINING_PROTEINS, TESTING_PROTEINS = getTrainTestByProteins(proteins_by_families, uniprot2pdb_mapping_dict, test_percentage=0.33)

        # Save protein family mapping
        with open(os.path.join(DATA_FOLDER, "proteins_by_families.pkl"), 'wb') as f:
            pickle.dump(proteins_by_families, f)
            
        with open(os.path.join(DATA_FOLDER, "protein_family_mapping.pkl"), 'wb') as f:
            pickle.dump(protein_family_mapping, f)
            
        with open(os.path.join(DATA_FOLDER, "uniprot2pdb_mapping_dict.pkl"), 'wb') as f:
            pickle.dump(uniprot2pdb_mapping_dict, f)

        with open(os.path.join(DATA_FOLDER, "training_proteins.pkl"), 'wb') as f:
            pickle.dump(TRAINING_PROTEINS, f)
            
        with open(os.path.join(DATA_FOLDER, "testing_proteins.pkl"), 'wb') as f:
            pickle.dump(TESTING_PROTEINS, f)

    training_data = IUPHAR(DATA_FOLDER, TRAINING_PROTEINS, training=True, pdb_id_index=1)
    testing_data = IUPHAR(DATA_FOLDER, TESTING_PROTEINS, training=False, pdb_id_index=1)
    
    # Load data and save keys
    if after_docking:
        load_and_save_data_by_keys(training_data.keys, training_data, training=True)
        load_and_save_data_by_keys(testing_data.keys, testing_data, training=False)
    
    else:
        if not os.path.exists("keys"):
            os.mkdir("keys")
        
        with open("keys/train_%s.pkl"%"IUPHAR", 'wb') as f:
            pickle.dump(training_data.keys, f)
            
        with open("keys/test_%s.pkl"%"IUPHAR", 'wb') as f:
            pickle.dump(testing_data.keys, f)
