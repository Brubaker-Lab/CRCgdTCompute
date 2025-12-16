This folder contains 4 scripts that are used for SCENIC-related analysis.

1. prepare_data_for_SCENIC.ipynb: Extract TIL and save AnnData as .loom.

2. SCENIC_hpc.py: On CWRU HPC, run this python script to run the pySCENIC workflow. The output of this step is adj.csv and reg.csv, and the latter will be used for further integration with NicheNet.

3. AfterSCENIC.ipynb: The jupyter notebook used for visualizing the SCENIC-predicted regulons specificity in each cell subtypes. This script also melt (like, pandas melt) the compacted SCENIC direct output of regulons into TF-target gene pairs in a .csv.

4. ChIPseq_SCENIC_verify_1.ipynb: The jupyter notebook used for getting the coordinates of +/- 5k bp from TSS of the target genes in the SCENIC prediction.

5. CHIPseq_SCENIC_verify_2.R: The R script to find overlaps of the public ChIP-seq peaks of TF around their SCENIC-predicted NicheNet undocumented target genes TSS, which is obtained after running ChIPseq_SCENIC_verify_1.ipynb