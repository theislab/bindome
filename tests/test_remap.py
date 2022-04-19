
import pytest
import bindome as bd

# move later

def test_remap_basic1(annpath):
    bd.constants.ANNOTATIONS_DIRECTORY = annpath
    df1 = bd.bindome.datasets.REMAP2020.get_remap_peaks('ZNF257')
    # assert len(ctf.datasets.example_peaks()) == 4
    assert df1.shape[0] > 100
    
def test_remap_basic2(annpath):
    bd.constants.ANNOTATIONS_DIRECTORY = annpath
    df1 = bd.bindome.datasets.REMAP2020.get_remap_peaks('FOXO1')
    # assert len(ctf.datasets.example_peaks()) == 4
    assert df1.shape[0] > 100