import os

import scanpy as sc


class scATAC:
    @staticmethod
    def PBMCs_10x_v2(datadir="../data/"):
        pass

        adata = sc.read_10x_h5(
            os.path.join(
                datadir,
                "10X/10k_Human_PBMCs_ATAC_v2_Chromium_Controller",
                "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5",
            ),
            gex_only=False,
        )
        adata.var["chr"] = adata.var_names.str.split(":").str[0]
        adata.var["start"] = adata.var_names.str.split(":").str[1].str.split("-").str[0].astype(int)
        adata.var["end"] = adata.var_names.str.split(":").str[1].str.split("-").str[1].astype(int)
        adata.var["peak.length"] = adata.var["end"] - adata.var["start"] + 1

        # reads = adata.X.data
        # fragments = (adata.X / 2).ceil()
        # get fragments
        adata.X = (adata.X / 2).ceil()

        return adata
