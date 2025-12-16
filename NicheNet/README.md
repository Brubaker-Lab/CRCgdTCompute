This folder contains 5 scripts that are used for NicheNet-related analysis. Please check them in the following order.

1. Check_SCENIC_GRN_source.R: To be honest, maybe this script is better to be in a NicheNet-joint-SCENIC folder. After getting the regulons result from SCENIC and formatting it using SCENIC>AfterSCENIC.ipynb, use this script to compare which pairs of SCENIC-identified TF-target genes are in NicheNet's gr_network_human_21122021.rds, which can be found in the publisher's repository.

2. Nichenet_SCENIC_integration.R: The R script to add the SCENIC-identified TF-target gene pairs to the NicheNet's default network and recompute the ligand-target matrix.

3. Result_NicheNet_DEG.R: The R script to run NicheNet ligand inference using the DEG (approach #1).

4. Result_Nichenet_geneformer.R: The R script to run NicheNet ligand inference using the Geneformer-identified perturbation targets (approach #2).

5. NicheNet_visualize.ipynb: The Python script to visualize the top ligands influencing the transition and calculate the enrichment of top ligands in CRC TME cell types.
