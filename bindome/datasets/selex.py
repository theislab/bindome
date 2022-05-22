
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
            if accession is not None and accession_id != accession:
                continue
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
    def load_tf_and_zero_reads(tf, data, **kwargs):
        data_sel_tf = data[(data['tf.name'] == tf)]  # & (grp['cycle'].astype(str) == '1')]
        data_sel_zero = data[(data['cycle'] == 0) & data['library'].isin(set(data[data['tf.name'] == tf][
                                                                              'library']))]  # & grp['accession'].isin(set(grp[grp['tf.name'] == tf]['accession']))]

        print(data_sel_tf.shape[0], data_sel_zero.shape[0])
        if data_sel_tf.shape[0] == 0 or data_sel_zero.shape[0] == 0:
            assert False
        # print('loading', tf, ':', library)
        reads_tf_next = SELEX.load_read_counts(tf, data=data_sel_tf, **kwargs)
        reads_zero_next = SELEX.load_read_counts(data=data_sel_zero, **kwargs)

        return reads_tf_next, reads_zero_next
    

    @staticmethod
    def load_read_counts(tf_name=None, data=None, library=None, fastq=True, k_skip=None):
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
            
            if k_skip is not None and k in k_skip:
                df_by_filename[k] = None
                continue
            
            # print(r['filename'], end='')
            p = os.path.join(selex_dir, r['accession'], r['filename'])

            i = 0
            seqlen = set()
            reads = []
            if fastq:
                for r in gzip.open(p):
                    if (i + 3) % 4 == 0:
                        # print(r)
                        s = r.strip().decode("utf-8")
                        if len(seqlen) == 0:
                            seqlen.add(len(s))
                        else:
                            if len(s) not in seqlen:
                                print(i, seqlen, len(s), s)
                        reads.append(s)
                    i += 1
                df = pd.DataFrame(reads, columns=['seq'])
            else:
                df = pd.read_csv(p)
                df.columns = ['seq']

            # force most frequent sequence is considered
            top_seqlen = df.seq.str.len().value_counts().sort_values().index[-1]
            # print(df.seq.str.len().value_counts().sort_values().index[-1])
            df = df[df['seq'].str.len() == top_seqlen]

            df = df.seq.value_counts().reset_index()
            df.columns = ['seq', 'counts']
            # print('# uniq reads/total counts %i/%i' % (df.shape[0], df['counts'].sum()))
            # df[df['tf.name'].str.contains('GATA')]

            # print(p, os.path.exists(p))
            uniq_seqlen = set(df['seq'].str.len())
            # print(uniq_seqlen)
            var_len_input = len(uniq_seqlen) != 1
            if var_len_input:
                for seqlen in uniq_seqlen:
                    print(seqlen)
                    print(df[df['seq'].str.len() == seqlen].head(1))
                assert var_len_input

            df_by_filename[k] = df
            # break
        return df_by_filename