import scanpy as sc
import decoupler as dc
import pandas as pd

def run_decoupler(ann_object_path, reg_path, out_path):
    ann_data = sc.read_h5ad(ann_object_path)
    reg = pd.read_csv(reg_path)

    dc.run_viper(mat=ann_data, net=reg, source='source', target='target', weight='weight', verbose=True, use_raw=False)

    estimates = ann_data.obsm['viper_estimate']
    estimates.to_csv(
        out_path + "decoupler_results.csv")
    print("Finished calculating TF activities")
