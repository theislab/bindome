import bindome as bd
import pandas as pd
import gzip
import os

class PBM():
    @staticmethod
    def gcn4_dream_v11():
        probound_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/pbm/freduce'
        df = pd.read_csv(os.path.join(probound_dir, 'gcn4.dream.v11.clean.txt'),
                        sep='\t')
        # print(df.head())
        df.columns = ['tf.name', 'tech', 'intensity', 'seq']
        return df