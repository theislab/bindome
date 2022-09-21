import numpy as np

import bindome as bd

# move later


def test_remap_basic1(annpath):
    bd.constants.ANNOTATIONS_DIRECTORY = annpath
    df1 = bd.bindome.datasets.REMAP2020.get_remap_peaks("ZNF257")
    # assert len(ctf.datasets.example_peaks()) == 4
    assert df1.shape[0] > 100


def test_remap_basic2(annpath):
    bd.constants.ANNOTATIONS_DIRECTORY = annpath
    df1 = bd.bindome.datasets.REMAP2020.get_remap_peaks("FOXO1")
    # assert len(ctf.datasets.example_peaks()) == 4
    assert df1.shape[0] > 100


def test_denormalize_bigwig(tmp_path):
    # TODO: Check that sum of counts is right
    from textwrap import dedent

    import pandas as pd

    path = tmp_path / "partial_SRX018625.tsv"
    with path.open("w") as f:
        # Writing part of SRX018625.bw as bed
        f.write(
            dedent(
                """
            Chromosome\tStart\tEnd\tValue
            chr1\t10005\t10016\t0.1185699999332428
            chr1\t10016\t10026\t0.2371390014886856
            chr1\t10026\t10028\t0.3557089865207672
            chr1\t10028\t10041\t0.592848002910614
            chr1\t10041\t10047\t0.47427898645401
            chr1\t10047\t10062\t0.592848002910614
            chr1\t10062\t10064\t0.47427898645401
            chr1\t10064\t10073\t0.3557089865207672
            chr1\t10073\t10083\t0.47427898645401
            chr1\t10083\t10100\t0.3557089865207672
            chr1\t10100\t10111\t0.2371390014886856
            chr1\t10111\t10119\t0.3557089865207672
            chr1\t10119\t10123\t0.47427898645401
            chr1\t10123\t10124\t0.592848002910614
            chr1\t10124\t10134\t0.47427898645401
            chr1\t10134\t10145\t0.592848002910614
            chr1\t10145\t10147\t0.47427898645401
            chr1\t10147\t10149\t0.3557089865207672
            chr1\t10149\t10151\t0.592848002910614
            chr1\t10151\t10152\t0.7114179730415344
            chr1\t10152\t10155\t0.8299869894981384
            chr1\t10155\t10156\t0.7114179730415344
            chr1\t10156\t10157\t0.8299869894981384
            chr1\t10157\t10159\t0.9485570192337036
            chr1\t10159\t10160\t0.8299869894981384
            chr1\t10160\t10170\t0.9485570192337036
            chr1\t10170\t10176\t0.8299869894981384
            chr1\t10176\t10185\t0.9485570192337036
            chr1\t10185\t10187\t0.7114179730415344
            chr1\t10187\t10188\t0.592848002910614
            chr1\t10188\t10192\t0.47427898645401
            chr1\t10192\t10193\t0.3557089865207672
            chr1\t10193\t10196\t0.2371390014886856
            chr1\t10196\t10248\t0.1185699999332428
            chr1\t10356\t10392\t0.1185699999332428
            chr1\t10393\t10406\t0.1185699999332428
            chr1\t10406\t10410\t0.2371390014886856
            chr1\t10410\t10429\t0.3557089865207672
            chr1\t10429\t10430\t0.2371390014886856
            """
            )
        )
    df = pd.read_table(str(path))

    norm_factor = bd.datasets.ChIPAtlas.normalization_factor("SRX018625", "hg38")
    denormalized = df.Value / norm_factor

    # Check that outputs are roughly integers
    rounded = np.around(denormalized.unique(), decimals=2)
    expected = np.arange(1, 10, dtype=float)

    assert np.all(np.isin(rounded, expected))
