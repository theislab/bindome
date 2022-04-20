
import bindome as bd
import pandas as pd
import gzip
import os

class SELEX():
    
    @staticmethod
    def get_data(accession='PRJEB3289'):
        selex_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/selex'
        os.path.exists(selex_dir)    
        
        filenames = [f for f in os.listdir(os.path.join(selex_dir, accession))]
        tfs = set([f.split('_')[0].upper() for f in filenames])
        # len(tfs)

        print('# filenames', len(filenames))
        data = pd.DataFrame(filenames, columns=['filename'])
        data['library'] = data['filename'].str.split('_').str[1 if accession == 'PRJEB3289' else 2]
        data['batch'] = data['filename'].str.split('_').str[2 if accession == 'PRJEB3289' else 1]
        data['cycle'] = data['filename'].str.split('_').str[3].str.split('.').str[0] # .astype(int)
        data['tf.name'] = [f.split('_')[0].upper() for f in filenames]
        # data.head()
        return data

    @staticmethod
    def load_read_counts(tf_name, data=None, library=None, accession='PRJEB3289'):
        if data is None:
            data = SELEX.get_data(accession=accession)

        # print(data.shape)
        # print(data.head())
        # print(data[data['tf.name'].str.contains(tf_name)])
        selex_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/selex'
        df_by_filename = {}
        
        data_sel = data[data['tf.name'].str.contains(tf_name)]
        # print(data_sel.shape)
        for ri, r in data_sel.iterrows():
            
            # print(r['library'])
            if library is not None and r['library'] != library:
                continue
            
            k = r['filename'].replace('.fastq.gz', '').replace('.txt.gz', '')
            # print(r['filename'], end='')
            p = os.path.join(selex_dir, accession, r['filename'])

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