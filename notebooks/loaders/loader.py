import yaml
from yaml.loader import SafeLoader
import anndata as ad
import numpy as np
import pandas as pd
import os
import gzip
import urllib
import bindome as bd

import warnings
warnings.filterwarnings('ignore')

class Loader():
    def __init__(self, data_dir, n_sample=None):
        self.data_dir = data_dir
        self.n_sample = n_sample
        self.filename = None
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)

    def load(self, url):
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)

        self.filename = url.split('/')[-1]
        if self.filename not in os.listdir(self.data_dir):
            print('Downloading from', url)
            urllib.request.urlretrieve(f'ftp://{url}', f'{self.data_dir}/{self.filename}')
            print('Done.')


class FASTQLoader(Loader):
    def __init__(self, data_dir, n_sample=None):
        super().__init__(data_dir, n_sample)

    def load(self, url):
        super().load(url)

        seqlen = set()
        reads = []

        # get every fourth line starting with index 1 (every second line in four line group)
        for i, r in enumerate(gzip.open(f'{self.data_dir}/{self.filename}')):
            if i  % 4 == 1:
                s = r.strip().decode("utf-8")
                seqlen.add(len(s))
                reads.append(s)
                if self.n_sample is not None and len(reads) >= self.n_sample:
                    break
        df = pd.DataFrame(reads, columns=["seq"])

        # only get those with the longest length
        df = df[df["seq"].str.len() == max(seqlen)]

        # occurance counts of each sequence
        df = df.seq.value_counts().reset_index()
        df.columns = ["seq", "counts"]

        # load into anndata
        adata = ad.AnnData(df['counts'].to_numpy().reshape(-1, 1)) 
        adata.obs_names = df['seq'].to_numpy()
        adata.var_names = ['count']

        return adata
    

    def read_yaml(self, yaml_path):
        with open(yaml_path) as f:
            data = yaml.load(f, Loader=SafeLoader)
            return data

    
    def fetch_ENA_urls(self, yaml_path):
        info = self.read_yaml(yaml_path)
        acc_ids = info['acc_ids']
        urls = []
        for id in acc_ids:
            if f'{id}.tsv' not in os.listdir('.'):
                print(f'Downloading {id}.tsv...')
                urllib.request.urlretrieve(id.join(info['download_meta_url']),
                                            filename=f'{self.data_dir}/{id}.tsv')

            with open(f'{self.data_dir}/{id}.tsv', 'r') as f:
                meta_df = pd.read_csv(f, delimiter='\t')
                # urls.append(meta_df['submitted_ftp'])
                urls += meta_df['submitted_ftp'].to_list()

        return urls