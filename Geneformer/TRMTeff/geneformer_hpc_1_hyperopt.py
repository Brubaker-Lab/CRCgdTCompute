#srun -p gpu --gres=gpu:1 --cpus-per-task=24 --mem=128G  --time=4200 --pty /bin/bash
import sys
import os
import pickle
from datetime import datetime
sys.path.append(os.getcwd())
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
from geneformer import TranscriptomeTokenizer
from geneformer import Classifier

storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/gdT_geneformer_TRMTeff'
output_prefix="gdT_activation"
vanilla_model = "/home/rxr456/Geneformer/gf-12L-95M-i4096"

tk = TranscriptomeTokenizer({"general_type": "general_type"}, nproc=39)
tk.tokenize_data(f"{storage_dir}", 
                 f"{storage_dir}",
                 "tokenized", 
                 file_format="h5ad")


cc = Classifier(classifier="cell",
                cell_state_dict = {"state_key": "general_type", "states": "all"},
                max_ncells=None,
                freeze_layers = 6,
                num_crossval_splits = 1,
                split_sizes = {"train": 0.6, "valid": 0.2, "test": 0.2},
                forward_batch_size=16,
                nproc=39)

cc.prepare_data(input_data_file=f"{storage_dir}/tokenized.dataset",
                output_directory=f"{storage_dir}/",
                output_prefix=output_prefix)

all_metrics = cc.validate(model_directory=vanilla_model,
                          prepared_input_data_file=f"{storage_dir}/{output_prefix}_labeled_train.dataset",
                          id_class_dict_file=f"{storage_dir}/{output_prefix}_id_class_dict.pkl",
                          output_directory=f"{storage_dir}/",
                          output_prefix=output_prefix,
                          n_hyperopt_trials=100,
                          predict_eval=True)