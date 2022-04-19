
import bindome as bd
import pandas as pd
import gzip
import os

class SELEX():
    
    @staticmethod
    def get_data():
        selex_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/selex'
        os.path.exists(selex_dir)    
        
        filenames = [f for f in os.listdir(os.path.join(selex_dir, 'PRJEB3289'))]
        tfs = set([f.split('_')[0].upper() for f in filenames])
        len(tfs)

        data = pd.DataFrame(filenames, columns=['filename'])
        data['library'] = data['filename'].str.split('_').str[1]
        data['batch'] = data['filename'].str.split('_').str[2]
        data['cycle'] = data['filename'].str.split('_').str[3].str.replace('.fastq.gz', '') # .astype(int)
        data['tf.name'] = [f.split('_')[0].upper() for f in filenames]
        # data.head()
        return data

    @staticmethod
    def load_read_counts(tf_name, data=None, library=None):
        if data is None:
            data = SELEX.get_data()

        selex_dir = bd.constants.ANNOTATIONS_DIRECTORY + '/selex'
        df_by_filename = {}
        for ri, r in data[data['tf.name'].str.contains(tf_name)].iterrows():
            
            if library is not None and r['library'] != library:
                continue
            
            k = r['filename'].replace('.fastq.gz', '')
            print(r['filename'], end='')
            p = os.path.join(selex_dir, 'PRJEB3289', r['filename'])

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