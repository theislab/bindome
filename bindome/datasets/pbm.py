import bindome as bd
import pandas as pd
import gzip
import os
import anndata

class PBM():
    @staticmethod
    def gcn4_dream_v11():
        fredure_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/pbm/freduce'
        df = pd.read_csv(os.path.join(freduce_dir, 'gcn4.dream.v11.clean.txt'),
                        sep='\t')
        # print(df.head())
        df.columns = ['tf.name', 'tech', 'intensity', 'seq']
        return df

    @staticmethod
    def uniprobe():
        h5ad_path = os.path.join(bd.constants.ANNOTATIONS_DIRECTORY, 'pbm', 'uniprobe', 'All_deBruijn.h5ad')
        if not os.path.exists(h5ad_path):
            print(h5ad_path)
            assert os.path.exists(h5ad_path)
        ad = anndata.read_h5ad(h5ad_path)
        return ad
