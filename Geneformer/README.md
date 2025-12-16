This folder contains 2 scripts used outside geneformer and a folder of 7 scripts used with geneformer that are used for Geneformer-related analysis.

1. Geneformer_prepare.ipynb: convert Gene symbols to ENSEMBL gene IDs (ENSG), drop Geneformer-unrelated attribute in the AnnData object, and save a .h5ad object as the input of Geneformer.

2. TRMTeff>geneformer_hpc_1_hyperopt.py: The actual python script to run Geneformer functions to tokenize the AnnData object and finetune their pretrained model with hyperparameter optimization.

3. TRMTeff>geneformer_hpc_1_printemb.py: The actual python script to run Geneformer functions to plot the classification performance of the finetuned model.

4. TRMTeff>geneformer_hpc_2_state_emb.py: The actual python script to run Geneformer functions to perform in silico perturbation.

5. TRMTeff>geneformer_hpc_2_state_emb.slurm: The script to submit the job of running geneformer_hpc_2_state_emb.py to CWRU HPC

6. TRMTeff>geneformer_hpc_2_state_emb_debug.slurm: The script to test if geneformer_hpc_2_state_emb.py can run bug-free before sumbitting it to HPC, which is very hard to get into the queue because of the computational resources it needs.

7. TRMTeff>geneformer_slurm_1_hyperopt.slurm: The script to submit the job of running geneformer_hpc_1_hyperopt.py to CWRU HPC
 
8. TRMTeff>geneformer_slurm_1_printemb.slurm: The script to submit the job of running geneformer_hpc_1_printemb.py to CWRU HPC

9. Geneformer_visualization.ipynb: The python script to visualize the Geneformer results. To visualize which significant targets are also DEGs, please run DEG>DEG.ipynb from another folder first.