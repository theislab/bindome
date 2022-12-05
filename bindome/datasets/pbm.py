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
    def pbm_paralogs():
        pbm_data_dir = os.path.join(bd.constants.ANNOTATIONS_DIRECTORY, 'pbm', 'GSE97794')
        print(pbm_data_dir)
        df = []

        sample_info =  pd.read_csv(os.path.join(pbm_data_dir, 'description.tsv'), sep='\t')

        # print(sample_info.head())
        sample_info = sample_info.set_index('sample')['description'].to_dict()

        for f in os.listdir(pbm_data_dir):
            if f == 'description.tsv' or '8mers' in f:
                continue

            df2 = pd.read_csv(os.path.join(pbm_data_dir, f))
            # print(df2.shape)
            df2['filename'] = f
            df2['description'] = sample_info[f.split('_')[0]]
            df2.append(df)
            df.append(df2)

        df = pd.concat(df)
        return df

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
