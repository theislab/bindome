
import bindome as bd
import pandas as pd
import gzip
import os
import numpy as np

class SELEX():
    
    @staticmethod
    def get_data(accession=None):
        selex_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/selex'
        os.path.exists(selex_dir)            
        filenames = []
        tfs = set()
        data_df = []
        for accession_id in os.listdir(selex_dir):
            print(accession_id)
            filenames = [f for f in os.listdir(os.path.join(selex_dir, accession_id))]
            tfs = set([f.split('_')[0].upper() for f in filenames]).union(tfs)
            # len(tfs)
            print('# filenames', len(filenames))
            data = pd.DataFrame(filenames, columns=['filename'])
            
            opt_library = 1 if accession_id == 'PRJEB3289' else (-1 if accession_id == 'PRJEB9797' else -2 if accession_id == 'PRJEB14744' else 2)
            data['library'] = data['filename'].str.split('.').str[0].str.split('_').str[opt_library]
            opt_batch = 2 if accession_id == 'PRJEB3289' else (1 if accession_id == 'PRJEB9797' else 1 if accession_id == 'PRJEB14744' else 2)
            data['batch'] = data['filename'].str.split('.').str[0].str.split('_').str[opt_batch]
            opt_cycle = -1 if accession_id == 'PRJEB3289' else (2 if accession_id == 'PRJEB9797' else 1 if accession_id == 'PRJEB14744' else 2)
            data['cycle'] = np.where(data['filename'].str.contains('Zero'), 0,
                                     data['filename'].str.split('.').str[0].str.split('_').str[opt_cycle])
            data['tf.name'] = [f.split('_')[0].upper() for f in filenames]
            data['accession'] = accession_id
            data_df.append(data)
        # data.head()
        return pd.concat(data_df).reset_index(drop=True)
    
    

    @staticmethod
    def load_read_counts(tf_name=None, data=None, library=None):
        if data is None:
            data = SELEX.get_data()

        # print(data.shape)
        # print(data.head())
        # print(data[data['tf.name'].str.contains(tf_name)])
        selex_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/selex'
        df_by_filename = {}
        
        data_sel = data[data['tf.name'].str.contains(tf_name)] if tf_name is not None else data
        
        # print(data_sel.shape)
        for ri, r in data_sel.iterrows():
            
            # print(r['library'])
            if library is not None and r['library'] != library:
                continue
            
            k = r['filename'].replace('.fastq.gz', '').replace('.txt.gz', '')
            # print(r['filename'], end='')
            p = os.path.join(selex_dir, r['accession'], r['filename'])

            i = 0
            reads = []
            for r in gzip.open(p):
                if (i + 3) % 4 == 0:
                    # print(r)
                    reads.append(r.strip().decode("utf-8"))
                i += 1
            df = pd.DataFrame(reads, columns=['seq'])
            df = df.seq.value_counts().reset_index()
            df.columns = ['seq', 'counts']
            print('# uniq reads/total counts %i/%i' % (df.shape[0], df['counts'].sum()))
            # df[df['tf.name'].str.contains('GATA')]

            df_by_filename[k] = df
            # break
        return df_by_filename