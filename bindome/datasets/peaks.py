import pandas as pd
import pybedtools
from pybedtools import BedTool


def get_overlapping_peaks(df1, df2):
    a = df1.copy()
    a = a[["chr", "start", "end", "coordinate"]]
    a.columns = ["chr", "start", "end", "k"]

    b = df2.copy()
    b = b[["chr", "start", "end", "coordinate"]]
    b.columns = ["chr", "start", "end", "k"]

    bed_a = BedTool.from_dataframe(a)
    bed_b = BedTool.from_dataframe(b)
    ab = bed_a.intersect(bed_b, wao=True).to_dataframe()
    pybedtools.cleanup()  # remove_all=True)

    ab["has.intersect"] = ab["itemRgb"] > 0
    ab["has.intersect"].value_counts()

    # print(a['k'].isin(ab['thickEnd']).value_counts())
    # print(b['k'].isin(ab['thickEnd']).value_counts())

    a["is.merged"] = a["k"].isin(set(ab[ab["has.intersect"]]["name"]))
    b["is.merged"] = b["k"].isin(set(ab[ab["has.intersect"]]["thickEnd"]))

    a[~a["is.merged"]]
    b[~b["is.merged"]]

    overlap = ab[ab["has.intersect"]]
    common1 = df1[df1["coordinate"].isin(overlap["name"])]
    df2[df2["coordinate"].isin(overlap["thickEnd"])]
    c = BedTool.from_dataframe(
        pd.concat(
            [common1[["chr", "start", "end", "coordinate"]], common1[["chr", "start", "end", "coordinate"]]]
        ).sort_values(["chr", "start"], ascending=[True, True])
    )

    merged = None
    if c.merge(header=True).count() > 0:
        merged = c.merge(header=True).to_dataframe()
        merged.columns = ["chr", "start", "end"]
        merged["coordinate"] = merged["chr"] + ":" + merged["start"].astype(str) + "-" + merged["end"].astype(str)
        merged["3"] = "merged"
    else:
        merged = pd.DataFrame(columns=["chr", "start", "end", "3"])

    out = pd.concat([df1, df2, merged])
    out = out.rename(columns={"3": "label"})

    out = out[["chr", "start", "end", "coordinate", "label"]]
    return overlap, merged, out
