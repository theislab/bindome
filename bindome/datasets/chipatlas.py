import os

import pandas as pd

import bindome as bd


def _parse_record(element) -> dict:
    """Retrieve record info from element."""
    from html import unescape

    out = {}
    first_pass = dict(
        zip(
            [x.text for x in element.findAll("dt")],
            [unescape(x.text) for x in element.findAll("dd")],
        )
    )

    out["Number of total reads"] = int(first_pass["Number of total reads"])
    out["Reads aligned (%)"] = float(first_pass["Reads aligned (%)"])
    out["Duplicates removed (%)"] = float(first_pass["Duplicates removed (%)"])
    out["Number of peaks"] = first_pass["Number of peaks"]

    return out


def _get_record_elements(soup):
    return [x.parent.parent for x in soup.findAll(string="Number of total reads")]


def _get_sequencing_metadata_from_page(soup):
    metadata_element = soup.findAll(attrs={"class": "mdata"})[-1]
    records = [_parse_record(x) for x in _get_record_elements(metadata_element)]
    keys = [x.text for x in metadata_element.findAll(attrs={"class": "small"})]
    return dict(zip(keys, records))


def _get_sequencing_metadata(chip_atlas_id: str) -> dict:
    import requests
    from bs4 import BeautifulSoup

    url = f"https://chip-atlas.org/view?id={chip_atlas_id}"
    req = requests.get(url, "html.parser")
    return _get_sequencing_metadata_from_page(BeautifulSoup(req.text, "html"))


class ChIPAtlas:
    @staticmethod
    def get_db_path(**kwargs):
        return (
            os.path.join(bd.constants.ANNOTATIONS_DIRECTORY, "epigenomes/chipatlas")
            if kwargs.get("dbpath") is None
            else kwargs.get("dbpath")
        )

    @staticmethod
    def get_experiments_list(update_columns=False, dbpath=None, **kwargs):

        p = os.path.join(ChIPAtlas.get_db_path(dbpath=dbpath), "experimentList.tab.gz")
        if update_columns:
            print("update colummns=True. Please check whether this is necessary before removing flag.")
            # check longest line to add columns
            longest_line = [len(s.split("\t")) for s in open(p)]
            headers = list(range(1, max(longest_line) + 1))
            lines = ["\t".join(map(str, headers))] + [s for s in open(p)]
            writer = open(p, "w")
            for line in lines:
                writer.write(line)
            writer.close()

        df = pd.read_csv(p, sep="\t")
        schema = ChIPAtlas.get_experiments_list_schema(dbpath=dbpath)
        df.columns = list(schema["Description"][:9]) + ["metadata.%i" % i for i in range(1, df.shape[1] - 9 + 1)]
        return df

    @staticmethod
    def get_experiments_list_schema(dbpath=None):
        path = os.path.join(ChIPAtlas.get_db_path(dbpath=dbpath), "experimentList_schema.tab")
        print("reading", path)
        return pd.read_csv(path, sep="\t")

    @staticmethod
    def normalization_factor(exp_id, genome_version) -> float:
        metadata = _get_sequencing_metadata(exp_id)[genome_version]
        n_reads = metadata["Number of total reads"]
        pct_aligned = metadata["Reads aligned (%)"]
        metadata["Duplicates removed (%)"]

        # TODO looks way more like integers without removing dupes
        return 1 / (n_reads * (pct_aligned / 100) / 1_000_000)  # * (1 - pct_dupes_removed / 100)

    @staticmethod
    def get_target_genes_local(genome, tf_name, distance_kbp=1, output_dir=None, download=True):
        http_path = "http://dbarchive.biosciencedbc.jp/kyushu-u/%s/target/%s.%i.tsv" % (
            genome,
            tf_name,
            distance_kbp,
        )
        if output_dir is None:
            output_dir = os.path.join(bd.constants.ANNOTATIONS_DIRECTORY, "chipseq/chipatlas/target")

        if not os.path.exists(output_dir):
            print(os.path.exists(output_dir), output_dir)
            print("please setup output dir")
            return
        output_dir = os.path.join(output_dir, genome)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        output_path = os.path.join(output_dir, basename(http_path))
        print(os.path.exists(output_path), output_path)
        if not os.path.exists(output_path) and download:
            print("wget", http_path, output_path)

            system("wget " + http_path + " -O " + output_path)
            if filesize(output_path) == 0:
                remove(output_path)
        try:
            return DataFrameAnalyzer.read_tsv(output_path)
        except OSError:
            print("preblem reading file", IOError.message)
            return None

    @staticmethod
    def get_chip_atlas_by_species():
        chipatlas = ChIPAtlas.get_experiments_list()
        chipatlas["species"] = chipatlas["Genome assembly"].map(
            {
                "hg19": "human",
                "mm9": "mouse",
                "sacCer3": "yeast",
                "dm3": "fly",
                "ce10": "worm",
                "rn6": "zebrafish",
            }
        )
        chipatlas_by_sp = {sp: grp for sp, grp in chipatlas.groupby("species")}
        return chipatlas_by_sp

    @staticmethod
    def get_peaks(genome, experiment_id, datadir=None, peaks_thr=5):
        if peaks_thr not in {10, 20, 5}:
            print(peaks_thr, "not valid as a peaks_thr...")
            stop()

        bed_bkp_path = None
        datadir = os.path.join(ChIPAtlas.get_db_path(), "peaks", genome)
        if not os.path.exists(datadir):
            os.mkdir(datadir)
        bed_bkp_path = os.path.join(datadir, experiment_id + "_" + str(peaks_thr) + ".bed.gz")
        if os.path.exists(bed_bkp_path):
            print("loading from saved BED", bed_bkp_path)
            return pd.read_csv(bed_bkp_path, sep="\t")

        peaks_thr = str(peaks_thr).zfill(2)
        p = "http://dbarchive.biosciencedbc.jp/kyushu-u/{}/eachData/bed{}/{}.{}.bed".format(
            genome,
            peaks_thr,
            experiment_id,
            peaks_thr,
        )
        print("querying in ChIP-atlas (peak thr (-log10(Q)=%s)..." % peaks_thr)
        # print(p)
        bed = pd.read_csv(p, sep="\t", header=None)
        if bed.shape[0] == 0:
            print("empty dataframe")
            return None
        bed.columns = [
            "chr",
            "start",
            "end",
            "id",
            "score",
            ".",
            "fold.change",
            "minus.log10.pval",
            "minus.log10.qval",
            "summit",
        ]

        # print(bed.head())
        bed["k"] = bed["chr"] + ":" + bed["start"].astype(str) + "-" + bed["end"].astype(str)
        # from lib.SequenceMethods import SequenceMethods
        # bed = SequenceMethods.parse_range2coordinate(bed, ['chr', 'start', 'end'], 'k.summit')
        if bed_bkp_path is not None:
            print("saving to output...", bed_bkp_path)
            bed.to_csv(bed_bkp_path, sep="\t")
        return bed

    @staticmethod
    def get_target_genes(genome, tf_name, distance_kbp=1):
        if distance_kbp not in {1, 5, 10}:
            print(distance_kbp, "not valid as a distance threshold...")
            stop()

        p = "http://dbarchive.biosciencedbc.jp/kyushu-u/%s/target/%s.%i.tsv" % (
            genome,
            tf_name,
            distance_kbp,
        )
        print("querying target genes in ChIP-atlas (distance thr = %iKbp)..." % distance_kbp)
        print(p)
        df = pd.read_csv(p, sep="\t")
        return df
