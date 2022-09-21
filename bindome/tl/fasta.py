import gzip
import random
import tempfile
from os import system
from os.path import exists, join
from random import shuffle

import numpy as np
import pandas as pd


def trim_fasta_sequence(fasta_path, bp_width, output_path):
    """
    Take a fasta file and generate a new one that has trimmed sequences
    :param fasta_path:
    :param bp_width:
    :return:
    """
    output_rows = []
    fastas = get_fastas_from_file(fasta_path, uppercase=True)
    print(len(fastas))

    for header, seq in fastas:
        for i in range(len(seq) / bp_width):
            subsequence = seq[i * bp_width : (i + 1) * bp_width]
            if len(subsequence) == bp_width:
                output_rows.append([header, str(i), subsequence])

    # write protein to output
    writer = open(output_path, "w")
    for r in output_rows:
        writer.write(">" + r[0] + "_" + r[1] + "\n" + r[2] + "\n")
    writer.close()


def convert_bed_to_bed_max_position(bed_peaks_path, bed_peaks_output_path, compression=None):
    """
    Convert a BED6+4 into a BED file that indicates the position of the peak
    max only
    :param bed_peaks_output_path:
    :return:
    """
    # create a new coordinates file with flanking sequences
    df = pd.read_csv(
        bed_peaks_path,
        sep="\t",
        index_col=False,
        names=["chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"],
    )
    print("here...")
    df["startPeak"] = df["chromStart"] + df["peak"]
    df["endPeak"] = df["startPeak"] + 1
    df["id"] = df["chrom"].astype(str) + ":" + df["startPeak"].astype(str) + "-" + df["endPeak"].astype(str)
    df = df[["chrom", "startPeak", "endPeak", "id"]]
    print("saving tmp file...")
    df.to_csv(bed_peaks_output_path, header=False, sep="\t", index=False, compression=compression)


def convert_fastq_to_fasta(fastq_path, fasta_path):
    fastq_analyzer = FastQAnalyzer()
    seqs = fastq_analyzer.get_sequences(fastq_path)
    write_fasta_from_sequences(seqs, fasta_path)


def convert_fasta_to_bed(fasta_path, bed_path):
    headers = [fa[0] for fa in get_fastas_from_file(fasta_path)]
    writer = open(bed_path, "w")
    for h in headers:
        chromosome, peak_range = h.split(":")
        start, end = peak_range.split("-")
        writer.write("\t".join([chromosome, start, end]) + "\n")
    writer.close()


def convert_bed_to_peaks_from_summit(bed_path, bp_flanking=50, stop_at=None):
    """

    :param bed_path: The path to our BED file
    :param output_path: The output bed that will be created
    :param bp_flanking: If use_peak is True, then flanking regions will
    (See https://www.biostars.org/p/102710/ for format description
    be calculated from this file
    :return:
    """
    # create a new coordinates file with flanking sequences
    print("reading tmp bed file...")
    df = pd.read_csv(
        bed_path,
        sep="\t",
        index_col=False,
        names=["chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"],
        nrows=stop_at,
    )
    print("here...")
    df["startFromPeak"] = df["chromStart"] + df["peak"] - bp_flanking
    df["endFromPeak"] = df["chromStart"] + df["peak"] + bp_flanking
    df = df[["chrom", "startFromPeak", "endFromPeak"]]
    tmp_bed_path = tempfile.mkstemp()[1]
    print("saving tmp file...")
    df.to_csv(tmp_bed_path, header=False, sep="\t", index=False)
    return tmp_bed_path


def scrambled_fasta_order(p, tmp_path=None, random_seed=None):
    fasta = get_fastas_from_file(p)
    if random_seed is not None:
        random.seed(random_seed)
    random.shuffle(fasta)
    tmp_path = tempfile.mkstemp()[1] if tmp_path is None else tmp_path
    write_fasta_from_sequences(fasta, tmp_path, uppercase=True)
    return tmp_path


def randomize_fasta(p, tmp_path=None):
    fasta = get_fastas_from_file(p)

    seqs = [randomize_sequence(s[1]) for s in fasta]
    tmp_path = tempfile.mkstemp()[1] if tmp_path is None else tmp_path
    write_fasta_from_sequences(seqs, tmp_path)
    return tmp_path


def get_sequences_from_bed(bed_path_or_dataframe, genome="hg19", **kwargs):
    fasta_path = tempfile.mkstemp()[1] if kwargs.get("fasta_path") is None else kwargs.get("fasta_path")
    if "fasta_path" in kwargs:
        del kwargs["fasta_path"]
    convert_bed_to_fasta(bed_path_or_dataframe, fasta_path, genome=genome, **kwargs)
    return get_fastas(fasta_path, **kwargs)


def convert_bed_to_fasta(bed_path_or_dataframe, fasta_path, genome="hg19", **kwargs):

    bed_path = None
    overwrite = kwargs.get("overwrite", True)
    if overwrite or not exists(fasta_path):
        if isinstance(bed_path_or_dataframe, pd.DataFrame):
            bed_path_tmp = kwargs.get("bed_path", tempfile.mkstemp()[1])
            bed_path_or_dataframe[list(bed_path_or_dataframe.columns)[:3]].to_csv(
                bed_path_tmp, header=None, sep="\t", index=None
            )
            bed_path = bed_path_tmp
        else:
            bed_path = bed_path_or_dataframe

        print(bed_path_tmp)
        print("genome", genome, genome == "mm10")
        gen_path = kwargs.get("gen_path", None)
        if gen_path is None:  # a genome was not given as input
            print("options")
            if genome == "hg19" or genome == "hg38":
                gen_path = "/g/scb2/zaugg/zaugg_shared/annotations/hg19/referenceGenome/%s.fa" % genome
                if not exists(gen_path):
                    gen_path = f"C:\\Users\\ignacio\\Dropbox\\{genome}\\{genome}.fa"
                    if not exists(gen_path):
                        gen_path = "/mnt/c/Users/ignacio.ibarra/Dropbox/annotations/{}/genome/{}.fa".format(
                            genome, genome
                        )
                    assert exists(gen_path)
            elif genome == "mm10":
                gen_path = join(
                    kwargs.get("basedir", "/g/scb2/zaugg/zaugg_shared/annotations/mm10/reference"), "mm10.fa"
                )
            elif genome == "mm9":
                gen_path = "/g/scb2/zaugg/zaugg_shared/annotations/mm9/Bowtie2/mm9/mm9.fa"
        print(gen_path)
        convert_bed_to_fasta_hg19(bed_path, fasta_path, genome_path=gen_path, **kwargs)


def convert_bed_to_fasta_hg19(bed_path_or_df, fasta_path, genome_path=None, force_strand=False, **kwargs):
    """

    :param bed_path: The path to our BED file
    :param fasta_path: The output fasta that will be created
    :param use_peak_max: If True, we will extract w.r.t. to peak position (+/- 50bp default)
    (See https://www.biostars.org/p/102710/ for format description
    :param bp_flanking: If use_peak is True, then flanking regions will
    be calculated from this file
    :return:
    """
    # create a new coordinates file with flanking sequences

    bed_path = None
    if not isinstance(bed_path_or_df, str):
        tmp_bed = tempfile.mkstemp()[1]
        DataFrameAnalyzer.to_tsv(bed_path_or_df, tmp_bed, header=False)
        bed_path = tmp_bed
    else:
        bed_path = bed_path_or_df

    n_peaks = kwargs.get("n", None)
    tmp_bed = kwargs.get("tmp_bed", None)
    use_peak_max = kwargs.get("use_peak_max", False)
    bp_flanking = kwargs.get("bp_flanking", 50)
    if use_peak_max:
        print("reading tmp bed file...")
        df = pd.read_csv(
            bed_path,
            sep="\t",
            index_col=False,
            names=[
                "chrom",
                "chromStart",
                "chromEnd",
                "name",
                "score",
                "strand",
                "signalValue",
                "pValue",
                "qValue",
                "peak",
            ],
        )

        if n_peaks is not None:
            df = df.head(n_peaks)

        print("here...")
        print(df.head())
        summit_index = kwargs.get("summit_index", None)
        if summit_index is None:
            df["startFromPeak"] = df["chromStart"] - bp_flanking
            df["endFromPeak"] = df["chromStart"] + bp_flanking
        else:
            df["startFromPeak"] = df["chromStart"] + df[df.columns[summit_index]] - bp_flanking
            df["endFromPeak"] = df["chromStart"] + df[df.columns[summit_index]] + bp_flanking

        if sum([k is None for k in df["peak"]]) != 0:
            df["startFromPeak"] = df["startFromPeak"] + df["peak"]
            df["startFromPeak"] = df["endFromPeak"] + df["peak"]

        df = df[["chrom", "startFromPeak", "endFromPeak"]]

        tmp_bed_path = tempfile.mkstemp()[1] if tmp_bed is None else tmp_bed
        print("saving tmp file...")
        df.to_csv(tmp_bed_path, header=False, sep="\t", index=False)
        bed_path = tmp_bed_path
        print("tmp bed file created...")

    if genome_path is not None:
        selected_path = genome_path

    print(exists(selected_path), selected_path)
    assert exists(selected_path)
    # bin = '/g/software/bin/bedtools'
    bin = "/g/funcgen/bin/bedtools"
    if not exists(bin):
        bin = "bedtools"

    args = [bin, "getfasta", "-fi", selected_path, "-bed", bed_path, "-fo", fasta_path] + (
        ["-s"] if force_strand else []
    )
    if kwargs.get("compress", False) or fasta_path.endswith(".gz"):
        args = (
            [bin, "getfasta", "-fi", selected_path, "-bed", bed_path]
            + (["-s"] if force_strand else [])
            + ["| gzip >", fasta_path]
        )

    # print bed_path
    # print fasta_path
    # print " ".join(args)
    print("running bedtools...")
    cmd = " ".join(args)
    print(cmd)
    system(cmd)


def create_bed_file(coordinates_table, bed_output_path):
    """
    CHROMOSOME_ID    START    END
    """

    if isinstance(coordinates_table, pd.DataFrame):
        coordinates_table = [r.values[:3] for ri, r in coordinates_table.iterrows()]

    # write HOMER motifs to tmp BED file. This one will be used later for validation
    is_gzip = bed_output_path.endswith(".gz")
    writer = gzip.open(bed_output_path, "w") if is_gzip else open(bed_output_path, "w")
    for r in coordinates_table:
        writer.write("\t".join(map(str, [r[0], r[1], r[2]])) + "\n")
    writer.close()


def get_fastas(fasta_path, uppercase=False, stop_at=None, na_remove=False, is_gzip=False, **kwargs):
    return get_fastas_from_file(fasta_path, uppercase=uppercase, stop_at=stop_at, is_gzip=is_gzip, **kwargs)


def get_fastas_from_file(fasta_path, uppercase=False, stop_at=None, na_remove=False, is_gzip=False, **kwargs):
    fastas = []
    seq = None
    header = None
    for r in gzip.open(fasta_path, mode="rt") if (fasta_path.endswith(".gz") or is_gzip) else open(fasta_path):
        r = r.strip()
        if r.startswith(">"):
            if seq != None and header != None:
                fastas.append([header, seq])
                if stop_at != None and len(fastas) >= stop_at:
                    break
            seq = ""
            header = r[1:]
        else:
            if seq != None:
                seq += r.upper() if uppercase else r
            else:
                seq = r.upper() if uppercase else r

    # append last fasta read by method
    if stop_at != None and len(fastas) < stop_at:
        fastas.append([header, seq])
    elif stop_at == None:
        fastas.append([header, seq])

    if kwargs.get("as_dict", False):
        if na_remove:
            return {h: s for h, s in fastas if not "N" in s}
        return {h: s for h, s in fastas}

    if na_remove:
        return [t for t in fastas if not "N" in t[1]]
    return fastas


def find_all(substring, dna):
    subs = [dna[i : i + len(substring)] for i in range(0, len(dna))]
    return [ind for ind, ele in enumerate(subs) if ele == substring]


def write_fasta_from_sequences(sequences, fasta_path, add_headers=True, uppercase=False):
    # print '# of sequences to write', len(sequences)
    writer = open(fasta_path, "w") if not fasta_path.endswith(".gz") else gzip.open(fasta_path, "w")
    # print 'writing fasta seq by seq...'
    for i, entry in enumerate(sequences):
        snp_id, seq = entry if len(entry) == 2 else [i, entry]
        next_seq = seq if not uppercase else seq.upper()
        s = ((">" + str(snp_id) + "\n") if add_headers else "") + (next_seq + "\n")
        writer.write(s)
    writer.close()


def concatenate(fasta_paths, output_path=None):
    output_path = tempfile.mkstemp()[1] if output_path is None else output_path
    with open(output_path, "w") as writer:
        for i, p in enumerate(fasta_paths):
            print(i, p)
            for r in open(p):
                writer.write(r)
    return output_path


def randomize_sequence(sequence):
    nucleotides = [nt for nt in sequence]
    shuffle(nucleotides)
    return "".join(nucleotides)


def get_gene_tss(upstream, downstream=0, genome="hg19"):
    assert genome == "hg19"

    all_ids_path = "/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/all_genes_hg19.tsv.gz"
    df = DataFrameAnalyzer.read_tsv_gz(all_ids_path)

    chromosome_ids = list(map(str, list(range(1, 23)))) + ["X", "Y"]
    df = df[df["chromosome_name"].isin({i for i in chromosome_ids})]
    df["start"] = np.where(
        df["strand"] == 1, df["transcription_start_site"] - upstream, df["transcription_start_site"] - downstream
    )
    df["start"] = np.where(df["start"] < 0, 0, df["start"])
    df["end"] = np.where(
        df["strand"] == 1, df["transcription_start_site"] + downstream, df["transcription_start_site"] + upstream
    )
    df["chromosome_name"] = "chr" + df["chromosome_name"]
    df = SequenceMethods.parse_range2coordinate(df, ["chromosome_name", "start", "end"], "range")
    return df
