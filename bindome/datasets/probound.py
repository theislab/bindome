import os

import pandas as pd

import bindome as bd


class ProBound:
    @staticmethod
    def imr90_gr():
        probound_dir = bd.constants.ANNOTATIONS_DIRECTORY + "/chipseq/probound"
        df = pd.read_csv(
            os.path.join(probound_dir, "countTable.0.IMR90_GR_chip-seq_rep1.tsv.gz"), sep="\t", header=None
        )
        df.columns = ["seq", 0, 1]
        return df

    @staticmethod
    def GR_mult_conc():
        probound_dir = bd.constants.ANNOTATIONS_DIRECTORY + "/chipseq/probound"
        df2 = []
        for i, conc in enumerate(["30", "300", "3000"]):
            df = pd.read_csv(
                os.path.join(probound_dir, "countTable.%i.GR_%s.tsv.gz" % (i, conc)), sep="\t", header=None
            )
            df.columns = ["seq", 0, 1]
            df["batch"] = conc
            df2.append(df)
        return pd.concat(df2)

    @staticmethod
    def ctcf(
        left_flank="ACACTCTTTCCCTACACGACGCTCTTCCGATCTTGACGTC",
        right_flank="GACGTCAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
        flank_length=7,
    ):
        probound_dir = bd.constants.ANNOTATIONS_DIRECTORY + "/chipseq/probound"
        df = pd.read_csv(os.path.join(probound_dir, "countTable.0.CTCF_r3.tsv.gz"), sep="\t", header=None)
        df.columns = ["seq", 0, 1]
        df["seq"] = left_flank[-flank_length:] + df["seq"] + right_flank[:flank_length]
        return df
