import bindome as bd
import pandas as pd
import gzip
import os

class ProBound():
    @staticmethod
    def imr90_gr():
        probound_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/chipseq/probound'
        df = pd.read_csv(os.path.join(probound_dir, 'countTable.0.IMR90_GR_chip-seq_rep1.tsv.gz'),
                        sep='\t', header=None)
        df.columns = ['seq', 0, 1]
        return df