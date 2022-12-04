import os

import anndata
import pandas as pd
import scipy.io

import bindome as bd


class PBM:
    @staticmethod
    def gcn4_dream_v11():
        fredure_dir = bd.constants.ANNOTATIONS_DIRECTORY + "/pbm/freduce"
        df = pd.read_csv(os.path.join(fredure_dir, "gcn4.dream.v11.clean.txt"), sep="\t")
        # print(df.head())
        df.columns = ["tf_name", "tech", "intensity", "seq"]
        return df

    @staticmethod
    def uniprobe():
        h5ad_path = os.path.join(bd.constants.ANNOTATIONS_DIRECTORY, "pbm", "uniprobe", "All_deBruijn.h5ad")
        if not os.path.exists(h5ad_path):
            print(h5ad_path)
            assert os.path.exists(h5ad_path)
        ad = anndata.read_h5ad(h5ad_path)
        return ad

    @staticmethod
    def pbm_homeo_affreg():

        matlab_path = os.path.join(bd.constants.ANNOTATIONS_DIRECTORY, "pbm", "affreg", "PbmDataHom6_norm.mat")
        mat = scipy.io.loadmat(matlab_path)

        data = mat["PbmData"][0]
        symbols = data[0][0][0][0][0]
        seqs_dbd = data[0][0][0][0][1]
        data[0][0][0][0][2]
        signal = data[0][1].T
        seqs = []
        names = []
        for i, s in enumerate(seqs_dbd):
            name = symbols[i][0][0] + "_" + str(i)
            seq = str(s[0][0])
            seqs += [seq]
            names.append(name)

        df = pd.DataFrame()
        df["seq"] = seqs
        df["name"] = names

        return df, signal
