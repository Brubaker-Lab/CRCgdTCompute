import scanpy as sc
import os
import glob
import numpy as np
import pandas as pd
# try the cancer only here
f_loom_path_scenic = '/home/rxr456/SCENIC_output/gdT_cancer_only.loom'
f_tfs = '/home/rxr456/SCENIC_database/hs_hgnc_curated_tfs.txt'
f_motif_path = '/home/rxr456/SCENIC_database/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
f_db_glob = '/home/rxr456/SCENIC_database/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
f_db_names = ' '.join( glob.glob(f_db_glob) )
adata = sc.read_loom(f_loom_path_scenic)

command = f"pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj_cancer_only.csv --num_workers 24"
os.system(command)

command = f"pyscenic ctx adj_cancer_only.csv {f_db_names} --annotations_fname {f_motif_path} --expression_mtx_fname {f_loom_path_scenic} --output reg_cancer_only.csv --mask_dropouts --num_workers 24"
os.system(command)

nGenesDetectedPerCell = np.sum(adata.X>0, axis=1)
nGenesDetectedPerCell = pd.DataFrame(nGenesDetectedPerCell)
percentiles = nGenesDetectedPerCell.quantile([.01, .05, .10, .50, 1])
print(percentiles)
f_pyscenic_output = "/home/rxr456/pyscenic_output_cancer_only.loom"
command = f"pyscenic aucell {f_loom_path_scenic} reg_cancer_only.csv --output {f_pyscenic_output} --num_workers 24"
os.system(command)