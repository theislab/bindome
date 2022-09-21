import os
from os import listdir
from os.path import exists, join

import pandas as pd

import bindome as bd


def get_parent_path(p):
    return os.path.abspath(os.path.join(os.path.dirname(p), ".."))


class REMAP2020:
    @staticmethod
    def get_remap_peaks(tf, genome="hg19", check=True, summit_extend=100, remove_non_canonical_chromosomes=True):

        if genome != "hg19":
            assert genome == "hg19"
        else:
            p = None
            if tf is not None:
                p = REMAP2020.get_remap_peaks_path(tf)
            else:
                p = "/g/scb2/zaugg/rio/data/remap2/ReMap2_nrPeaks_hg19.bed.gz"
            if p is None:
                print(tf, "path cannot be found. Please check name or REMAP database")
                return None

            # print('reading peak data from', p, exists(p))
            df = pd.read_csv(p, header=0 if "remap2" in p else None, sep="\t")

            # print('renaming columns')
            df.columns = ["chr", "start", "end"] + list(df.columns)[3:]

            # print('parsing coordinates')

            df["coordinate"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
            df["summit.start"] = df["6" if "6" in df else 6] - summit_extend
            df["summit.end"] = df["6" if "6" in df else 6] + summit_extend

            df = df[df["summit.start"] >= 0]

            if remove_non_canonical_chromosomes:
                # print('removing non canonical chromosomes')
                df = df[~df["chr"].str.contains("_")].reset_index(drop=True)

            df["k.summit"] = df["chr"] + ":" + df["summit.start"].astype(str) + "-" + df["summit.end"].astype(str)
            return df

    @staticmethod
    def separate_peaks_into_files(option="nrPeaks"):
        basedir = "/g/scb2/zaugg/rio/data/remap2"
        p = join(basedir, "ReMap2_$option_hg19.bed.gz").replace("$option", option)

        df = DataFrameAnalyzer.read_tsv_gz(p, header=None)
        output_dir = join(basedir, option)
        if not exists(output_dir):
            mkdir(output_dir)
        for motif_id, grp in df.groupby(3):
            output_path = join(output_dir, motif_id + ".tsv.gz")
            print(motif_id, grp.shape, output_path)
            DataFrameAnalyzer.to_tsv_gz(grp, join(output_dir, motif_id + ".tsv.gz"))

    @staticmethod
    def intersect_peaks_remap(tf1, tf2):
        df1 = REMAP2020.get_remap_peaks(tf1)
        df2 = REMAP2020.get_remap_peaks(tf2)

        df1 = df1[df1["chr"] != "chrM"].reset_index(drop=True)
        df2 = df2[df2["chr"] != "chrM"].reset_index(drop=True)

        out1 = PeaksAnalyzer.intersect_coordinates(
            df1[["chr", "summit.start", "summit.end", "k.summit"]],
            df2[["chr", "summit.start", "summit.end", "k.summit"]],
        )
        df3 = DataFrameAnalyzer.read_tsv(
            out1[-1], header=None, columns=["chr", "summit.start", "summit.end", "k.summit"]
        )
        df3 = SequenceMethods.parse_range2coordinate(df3, ["chr", "summit.start", "summit.end"], "k")
        out2 = PeaksAnalyzer.intersect_coordinates(
            df2[["chr", "summit.start", "summit.end", "k.summit"]],
            df1[["chr", "summit.start", "summit.end", "k.summit"]],
        )
        df4 = DataFrameAnalyzer.read_tsv(
            out2[-1], header=None, columns=["chr", "summit.start", "summit.end", "k.summit"]
        )
        df4 = SequenceMethods.parse_range2coordinate(df4, ["chr", "summit.start", "summit.end"], "k")

        return df1, df2, df3, df4

    @staticmethod
    def get_great_regions(tf_name):
        path = join(
            "/g/scb2/zaugg/rio/EclipseProjects/zaugglab/comb-TF-binding",
            "combinatorial_binding_analysis_invivo/data/GREAT_results_peak_gene_associations",
            "%s.tsv.gz" % tf_name,
        )
        assert exists(path)
        res = DataFrameAnalyzer.read_tsv_gz(path)
        res = SequenceMethods.parse_range2coordinate(res, ["seqnames", "start", "end"], "k")
        res = res.rename(columns={"k": "coordinate", "gene": "gene.name", "distTSS": "distance"})
        return res[["coordinate", "gene.name", "distance"]]

    @staticmethod
    def get_narrowpeaks_fasta(code, n_peaks, peak_size=200):
        narrowPeak_dir = "/g/scb2/zaugg/rio/data/remap2/narrowPeak"

        query_f = None
        for f in listdir(narrowPeak_dir):
            if code in f and f.endswith(".narrowPeak.gz"):
                query_f = f
                break
        assert query_f is not None

        p = join(narrowPeak_dir, f)
        hg19_peaks = p.replace(".narrowPeak.gz", ".bed")
        hg19_fa_path = hg19_peaks.replace(".bed", ".fa")
        if True or not exists(hg19_fa_path):
            df = DataFrameAnalyzer.read_tsv_gz(p, header=None)

            # sort by peak relevance
            df = df.sort_values(list(df.columns)[-2], ascending=False).head(min(df.shape[0], n_peaks))
            df["start.summit"] = df[1] + df[9] - peak_size / 2
            df["end.summit"] = df[1] + df[9] + peak_size / 2

            bed_tmp = tempfile.mkstemp()[1]
            DataFrameAnalyzer.to_tsv(df[[0, "start.summit", "end.summit"]], bed_tmp, header=None)

            # convert hg38 to hg19
            chain_hg38_to_hg19 = (
                "/g/scb2/zaugg/rio/EclipseProjects/zaugglab/invitro_vs_multiomics/data/hg38ToHg19.over.chain.gz"
            )
            cmd = [
                "/g/software/bin/liftOver",
                bed_tmp,
                chain_hg38_to_hg19,
                hg19_peaks,
                join("../../../data/unMapped_" + basename(hg19_peaks)),
            ]

            if not exists(hg19_peaks):
                system(" ".join(cmd))

            # last step is converting this bed file to a proper fasta
            if not exists(hg19_fa_path):
                FastaAnalyzer.convert_bed_to_fasta(hg19_peaks, hg19_fa_path, genome="hg19")

        return hg19_fa_path

    @staticmethod
    def get_tfs_by_dbd(dbd=None):
        tfs = REMAP2020.get_remap_ids()
        tf_data_bkp = "/g/scb2/zaugg/rio/data/remap2/tf_data.pkl"
        if not exists(tf_data_bkp):
            ensembl_by_gene = MyGeneAnalyzer.get_ensembl_by_symbol(tfs)
            DataFrameAnalyzer.to_pickle(ensembl_by_gene, tf_data_bkp)

        ensembl_by_gene = DataFrameAnalyzer.read_pickle(tf_data_bkp)
        dbd_by_ensembl = HumanTFs.get_dbd_by_ensembl()

        res = []
        for g in tfs:
            if g in ensembl_by_gene and ensembl_by_gene[g] in dbd_by_ensembl:
                res.append([g, ensembl_by_gene[g], dbd_by_ensembl[ensembl_by_gene[g]]])
            else:
                res.append([g, None, None])
        res = pd.DataFrame(res, columns=["gene.name", "ensembl", "DBD"])
        if dbd is not None:
            return res[res["DBD"] == dbd].reset_index(drop=True)
        return res

    @staticmethod
    def get_intersections_vs_remap(df, bkp_path=None, overwrite=False):

        print("loading ReMap2018 intersections...")
        if overwrite or bkp_path is None or not exists(bkp_path):
            ranges = None
            if len(df.columns) != 3:  # probably ooordinates: Convert into ranges
                ranges = SequenceMethods.parse_coordinate2range(df, df[df.columns[0]])
                ranges = ranges[ranges.columns[1:]]
            else:
                ranges = df
            print(ranges.head())

            peak_analyzer = PeaksAnalyzer()
            bed_path = tempfile.mkstemp()[1]

            if bkp_path is None:
                intersections_path = tempfile.mkstemp()[1]
            else:
                intersections_path = bkp_path
            DataFrameAnalyzer.to_tsv(ranges, bed_path, header=None)
            remap2_nrPeaks = "/g/scb2/zaugg/rio/data/remap2/ReMap2_nrPeaks_hg19.bed.gz"
            peak_analyzer.insersect_chip_seq_peaks(bed_path, remap2_nrPeaks, intersections_path)

        intersections_df = DataFrameAnalyzer.read_tsv(bkp_path, header=None)
        intersections_df.columns = [
            "chr",
            "start",
            "end",
            "chr.remap",
            "start.remap",
            "end.remap",
            "tf.remap",
            ".",
            "..",
            "summit.start",
            "summit.end",
            "shape",
            "overlap",
        ]
        intersections_df = intersections_df[intersections_df["overlap"] != 0]
        return intersections_df

    @staticmethod
    def get_remap_ids():
        return [
            tf.replace(".tsv.gz", "").replace(".bed.gz", "").upper()
            for tf in listdir(join(bd.constants.ANNOTATIONS_DIRECTORY, "remap/remap3/nrPeaks"))
            if tf.endswith(".bed.gz")
        ]

    @staticmethod
    def get_remap_peaks_path(tf_name):
        # tf_id = tf_name.lower()

        basedir = join(bd.constants.ANNOTATIONS_DIRECTORY, "remap/remap3/nrPeaks")

        if not exists(basedir):
            print("please include this data directory to run loading of peaks successfully")
            # print(exists(basedir), basedir)
            assert not exists(basedir)

        if "remap3" in basedir:
            p = join(basedir, tf_name + ".bed.gz")
        if "remap2" in basedir:
            p = join(basedir, tf_name + ".tsv.gz")

        # print(exists(p), p)
        if not exists(p):
            if not exists(p):
                # print tf_name, 'not found in REMAP2'
                return None
            assert exists(p)
        return p

    @staticmethod
    def write_peaks_as_bed(tf_name, output_basename):
        peaks = REMAP2020.get_remap_peaks(tf_name)
        DataFrameAnalyzer.to_tsv(
            peaks[["chr", "summit.start", "summit.end", "k.summit"]], output_basename + ".bed", header=None
        )

    @staticmethod
    def get_merged_peaks(extend=2500, remove_non_canonical_chromosomes=True):
        nr_peaks_path = "/g/scb2/zaugg/rio/data/remap2/ReMap2_nrPeaks_hg19.bed.gz"
        merged_peaks_path = nr_peaks_path.replace(".bed.gz", "_merged_overlap_%i.bed.gz" % (extend))
        if not exists(merged_peaks_path):
            cmd = "bedtools merge -d %i -i %s | gzip > %s" % (extend, nr_peaks_path, merged_peaks_path)
            print(cmd)
            system(cmd)

        # print('reading peak data from', merged_peaks_path)
        df = DataFrameAnalyzer.read_tsv_gz(merged_peaks_path, header=None)
        # print('renaming columns')
        df.columns = ["chr", "start", "end"] + list(df.columns)[3:]

        # print('parsing coordinates')
        print(df.head())
        df = SequenceMethods.parse_range2coordinate(df, ["chr", "start", "end"], "coordinate")
        df["summit.start"] = ((df["start"] + df["end"]) / 2 - extend).astype(int)
        df["summit.end"] = (df["summit.start"] + 2 * extend).astype(int)

        df = df[df["summit.start"] >= 0]

        if remove_non_canonical_chromosomes:
            # print('removing non canonical chromosomes')
            df = df[~df["chr"].str.contains("_")].reset_index(drop=True)

        df = SequenceMethods.parse_range2coordinate(df, ["chr", "summit.start", "summit.end"], "k.summit")
        return df

    @staticmethod
    def get_all_peaks(extend=100):
        return REMAP2020.get_remap_peaks(None, summit_extend=extend)

    @staticmethod
    def get_random_remap_peaks(n, ngroup=100, seed=10, query_tfs=None):
        tf_names = REMAP2020.get_remap_ids()
        if query_tfs is not None:
            tf_names = [tf for tf in tf_names if tf in query_tfs]
        random.seed(seed)
        final = []
        while sum([next.shape[0] for next in final]) < n:
            next_tf = tf_names[random.randint(0, len(tf_names) - 1)]
            if query_tfs is not None and not next_tf in query_tfs:
                continue
            next_df = REMAP2020.get_remap_peaks(next_tf)
            sample = next_df.sample(min(ngroup, next_df.shape[0]))
            print(next_tf, sample.shape[0])
            final.append(sample)
        final = pd.concat(final).reset_index(drop=True).sort_values(["chr", "start"], ascending=[True, True])
        return final

    @staticmethod
    def get_closest_tss_gene(coordinates, custom_coordinates=None, genome="hg19"):

        df = None
        if custom_coordinates is None and genome == "hg19":
            bed_tss_hg19 = None
            if genome == "hg19":
                bed_tss_hg19 = "/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/hg19_refseq_n_ensembl_n_pos.tsv"

            df = DataFrameAnalyzer.read_tsv(bed_tss_hg19)
            df["start"] = df["tss"].astype(int)
            df["end"] = (df["start"] + 1).astype(int)
            df = df.sort_values(["chr", "start"], ascending=[True, True])
        else:
            df = custom_coordinates

        out = PeaksAnalyzer.intersect_coordinates(
            coordinates[["chr", "start", "end"]],
            df[["chr", "start", "end", "SYMBOL", "ENSEMBL"]],
            mode="closest",
            additional_args=["-d"],
        )
        res = DataFrameAnalyzer.read_tsv(
            out[-1],
            header=None,
            columns=["chr", "start", "end", "chr.gene", "tss.start", "tss.end", "symbol", "ensembl", "distance"],
        )
        res = SequenceMethods.parse_range2coordinate(res, ["chr", "start", "end"], "k")
        return res
