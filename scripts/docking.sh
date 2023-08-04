export CUDA_VISIBLE_DEVICES=3

export HOME=esm/model_weights
export PYTHONPATH=$PYTHONPATH:/dfs/user/sttruong/DucWorkspace/DiffDock/esm

python datasets/esm_embedding_preparation.py \
--protein_ligand_csv ../graph_regression/iuphar/input_protein_ligand_True.csv \
--out_file data/prepared_for_esm_train.fasta


python esm/scripts/extract.py esm2_t33_650M_UR50D data/prepared_for_esm_train.fasta_0_None data/esm2_output \
--repr_layers 33 --include per_tok --truncation_seq_length 30000


python -m inference \
--protein_ligand_csv ../graph_regression/iuphar/input_protein_ligand_True.csv \
--out_dir results/training \
--inference_steps 20 \
--samples_per_complex 40 \
--batch_size 8 \
--start 3500 \
--end 3500

0-300: ampere3 - 1
300-600: ampere3 - 4
600-900: ampere3 - 3
900-1200: ampere3 - 7

1200-2700: ampere4 - 6
2700-3500: ampere3 - 1
3500-: ampere3 - 3


python datasets/esm_embedding_preparation.py \
--protein_ligand_csv ../graph_regression/iuphar/input_protein_ligand_False.csv \
--out_file data/prepared_for_esm_test.fasta


python esm/scripts/extract.py esm2_t33_650M_UR50D data/prepared_for_esm_test.fasta_0_None data/esm2_output \
--repr_layers 33 --include per_tok --truncation_seq_length 30000


python -m inference \
--protein_ligand_csv ../graph_regression/iuphar/input_protein_ligand_False.csv \
--out_dir results/testing \
--inference_steps 20 \
--samples_per_complex 40 \
--batch_size 8 \
--start 2700 \
--end 600

0 - 300: ampere4 - 0
300 - 600: ampere4 - 1
600 - 900: ampere4 - 4
900 - 1200: ampere4 - 5

1200 - 2700: ampere4 - 7
2700 - : ampere3 - 7
