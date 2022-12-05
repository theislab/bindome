import pandas as pd
import bindome as bd

def test_load_pbm_paralogs(annpath):
    bd.constants.ANNOTATIONS_DIRECTORY = annpath
    data = bd.bindome.datasets.PBM.pbm_paralogs()

    print(data['filename'].value_counts())
    print(data['description'].value_counts())
    assert isinstance(data, pd.DataFrame)
