import os

import numpy as np
import pandas as pd
import matplotlib.cbook as cbook

from .constants import BASES, AUC_STEP

class EvaluationReport():

    def __init__(self, df, modelname, output_path = None, overwrite = False):

        self.df = self.calculate_rates(df)
        self.output_path = output_path
        self.modelname = modelname
        self.overwrite = overwrite

    def calculate_rates(self, df):

        match_columns = list()
        mismatch_columns = list()
        insertion_columns = list()
        deletion_columns = list()
        for b1 in BASES:
            for b2 in BASES + ['-']:
                for b3 in BASES + ['-']:
                    for b4 in BASES:
                        if b2 == '-' and b3 == '-':
                            continue
                        if b2 == b3:
                            match_columns.append(b1 + b2 + '>' + b3 + b4)
                        else:
                            if b2 == '-':
                                deletion_columns.append(b1 + b2 + '>' + b3 + b4)
                            elif b3 == '-':
                                insertion_columns.append(b1 + b2 + '>' + b3 + b4)
                            else:
                                mismatch_columns.append(b1 + b2 + '>' + b3 + b4)

        self.signature_columns = match_columns+mismatch_columns+insertion_columns+deletion_columns

        # calculate the number of each event
        df.loc[:, 'matches'] = df[match_columns].sum(axis=1)
        df.loc[:, 'mismatches'] = df[mismatch_columns].sum(axis=1)
        df.loc[:, 'insertions'] = df[insertion_columns].sum(axis=1)
        df.loc[:, 'deletions'] = df[deletion_columns].sum(axis=1)
        # length of the alignment is equal to the sum of event plus two because the first
        # and last bases are not counted in these 3mer events
        df.loc[:, 'len_alignment'] = df[['matches', 'mismatches', 'insertions', 'deletions' ]].sum(axis=1)
        df.loc[:, 'len_alignment'] += 2

        # calculate rates by dividing counts by the length of the alignment
        df.loc[:, 'match_rate'] = df['matches']/df['len_alignment']
        df.loc[:, 'mismatch_rate'] = df['mismatches']/df['len_alignment']
        df.loc[:, 'insertion_rate'] = df['insertions']/df['len_alignment']
        df.loc[:, 'deletion_rate'] = df['deletions']/df['len_alignment']

        # annotate reads whose alignment is way too shorter relative to the length of the reference
        df.loc[df['len_alignment']/df['len_reference'] < 0.7, 'comment'] = 'shortalignment'

        # calculate homopolymer rates per base and in total
        df['total_homo_counts'] = df['homo_A_counts'] + df['homo_C_counts'] + df['homo_G_counts'] + df['homo_T_counts']
        df['total_homo_errors'] = df['homo_A_errors'] + df['homo_C_errors'] + df['homo_G_errors'] + df['homo_T_errors']
        df['total_homo_error_rate'] = df['total_homo_errors']/df['total_homo_counts']
        df['A_homo_error_rate'] = df['homo_A_errors']/df['homo_A_counts']
        df['C_homo_error_rate'] = df['homo_C_errors']/df['homo_C_counts']
        df['G_homo_error_rate'] = df['homo_G_errors']/df['homo_G_counts']
        df['T_homo_error_rate'] = df['homo_T_errors']/df['homo_T_counts']

        return df

    def to_csv(self, df, reportname):

        if self.output_path is None:
            return None
        else:
            ofile = os.path.join(self.output_path, 'report_' + self.modelname + '_' + reportname + '.csv')
            if os.path.isfile(ofile):
                if self.overwrite:
                    os.remove(ofile)
                else:
                    raise FileExistsError('Output file ({0}) already exists'.format(ofile))
            return df.to_csv(ofile, index=False, header = True)
                

    def calculate_boxplot_stats(self, df, columns):

        stats = ['mean', 'med', 'iqr', 'whishi', 'whislo', 'q1', 'q3']

        boxplotstats = {
            'model':list(),
            'stats':stats,
        }

        for col in columns:
            boxplotstats[col] = list()
            values = df[col]
            values = values[~np.isnan(df[col])]
            rate_stats = cbook.boxplot_stats(values)[0]
            for st in stats:
                boxplotstats[col].append(rate_stats[st])
                if col == columns[0]:
                    boxplotstats['model'].append(self.modelname)
        boxplotstats = pd.DataFrame(boxplotstats)
        return boxplotstats
        

    def count_abs_counts(self):

        reportname = 'absolutecounts'
        subdf = self.df[self.df['comment'] == 'pass']

        simple_report = {
            'model' : self.modelname,
            'total_bases_basecalled': int(subdf['len_basecalls'].sum()),
            'total_bases_reference': int(subdf['len_reference'].sum()),
            'total_bases_aligned': int(subdf['len_alignment'].sum()),
            'total_match_bases': int(subdf['matches'].sum()),
            'total_mismatch_bases': int(subdf['mismatches'].sum()),
            'total_insertion_bases': int(subdf['insertions'].sum()),
            'total_deletion_bases': int(subdf['deletions'].sum()),
            'total_homopolymer_bases_reference': int(subdf['total_homo_counts'].sum()),
            'total_homopolymer_bases_error': int(subdf['total_homo_errors'].sum()),
        }

        simple_report = pd.DataFrame(simple_report, index = [0])
        self.to_csv(simple_report, reportname)
        return simple_report

    def read_outcome_counts(self):

        reportname = 'readoutcomes'

        read_outcome_counts = dict(self.df['comment'].value_counts())
        read_outcome_counts['model'] = 'model1'
        read_outcome_counts = pd.DataFrame(read_outcome_counts, index = [0])
        cols = list(read_outcome_counts.columns)
        cols = [cols[-1]] + cols[:-1]
        read_outcome_counts = read_outcome_counts[cols]

        self.to_csv(read_outcome_counts, reportname)
        return read_outcome_counts

    def event_rates(self):

        reportname = 'eventrates'
        subdf = self.df[self.df['comment'] == 'pass']
        cols = [
            'match_rate', 
            'mismatch_rate', 
            'insertion_rate', 
            'deletion_rate'
        ]

        events = self.calculate_boxplot_stats(subdf, cols)
        self.to_csv(events, reportname)
        return events

    def homopolymer_rates(self):

        reportname = 'homopolymerrates'
        subdf = self.df[self.df['comment'] == 'pass']
        cols = [
            'total_homo_error_rate', 
            'A_homo_error_rate', 
            'C_homo_error_rate', 
            'G_homo_error_rate', 
            'T_homo_error_rate'
        ]

        events = self.calculate_boxplot_stats(subdf, cols)
        self.to_csv(events, reportname)
        return events

    def phredq_distributions(self):

        reportname = 'phredq'
        subdf = self.df[self.df['comment'] == 'pass']
        cols = [
            'phred_mean_correct', 
            'phred_mean_error', 
        ]

        events = self.calculate_boxplot_stats(subdf, cols)
        self.to_csv(events, reportname)
        return events

    def calculate_signatures(self):

        reportname = 'signatures'
        subdf = self.df[self.df['comment'] == 'pass']

        types = [
            "Match", 
            "Missmatch_A", 
            "Missmatch_C", 
            "Missmatch_G", 
            "Missmatch_T", 
            "Missmatch", 
            "Insertion", 
            "Deletion"
        ]

        signatures = list()
        for b in BASES:
            for b1 in BASES:
                for b2 in BASES:
                    for t in types:
                        signatures.append({
                            'Base': b, 
                            'Context': b1+b+b2, 
                            'Error': t, 
                            "Count": 0, 
                            "Rate": 0.0, 
                            "k": b+b1+b2+t
                        })
        signatures_df = pd.DataFrame(signatures)

        for s_col in self.signature_columns:
            k1 = s_col[1] if s_col[1] != '-' else s_col[3]
            k2 = s_col[0] + s_col[-1]
            if s_col[1] == '-':
                k3 = "Deletion"
            elif s_col[3] == '-':
                k3 = "Insertion"
            elif s_col[1] == s_col[3]:
                k3 = "Match"
            else:
                k3 = "Missmatch"
            p = np.where(signatures_df["k"] == k1+k2+k3)[0]
            signatures_df.iloc[p, 3] += np.sum(subdf[s_col])
            if k3 == 'Missmatch':
                k3 = "Missmatch_" + s_col[3]
                p = np.where(signatures_df["k"] == k1+k2+k3)[0]
                signatures_df.iloc[p, 3] += np.sum(subdf[s_col])

        for b in np.unique(signatures_df["Base"]):
            sub_df = signatures_df[signatures_df["Base"] == b]
            for c in np.unique(signatures_df["Context"]):
                subsub_df = sub_df[sub_df["Context"] == c]
                subsub_df = subsub_df[subsub_df["Error"].isin(types)]
                t = np.sum(subsub_df["Count"])
                signatures_df.loc[subsub_df.index, 'Rate'] = signatures_df.loc[subsub_df.index, 'Count']/t

        d_list = list()
        for b in np.unique(signatures_df["Base"]):
            sub_df = signatures_df[signatures_df["Base"] == b]
            sub_df = sub_df[sub_df["Error"].isin(types)]
            total_bases = np.sum(sub_df['Count'])
            for c in np.unique(signatures_df["Error"]):
                subsub_df = sub_df[sub_df["Error"] == c]
                total_errors = np.sum(subsub_df['Count'])
                d = {
                    'Base': 'Grouped', 
                    'Context': 'N'+b+'N', 
                    'Error': c, 
                    'Count': total_errors, 
                    'Rate': total_errors/total_bases, 
                    'k': b+c
                    }
                d_list.append(d)

        total_bases = np.sum(signatures_df['Count'])
        for c in np.unique(signatures_df["Error"]):
            sub_df = signatures_df[signatures_df["Error"] == c]
            total_errors = np.sum(sub_df['Count'])
            d = {
                'Base': 'Grouped', 
                'Context': 'NNN', 
                'Error': c, 
                'Count': total_errors, 
                'Rate': total_errors/total_bases, 
                'k': b+c
                }
            d_list.append(d)

        signatures_df = signatures_df.append(pd.DataFrame(d_list))
        signatures_df = signatures_df.drop(labels='k', axis = 1)
        signatures_df['model'] = self.modelname
        cols = list(signatures_df.columns)
        cols = [cols[-1]] + cols[:-1]
        signatures_df = signatures_df[cols]

        self.to_csv(signatures_df, reportname)

        return signatures_df

    def calculate_auc(self):
        
        reportname = 'auc'
        subdf = self.df[self.df['comment'] == 'pass']
        subdf = subdf.sort_values('phred_mean', ascending=True)
        subdf = subdf[~np.isnan(subdf['phred_mean'])]
        subdf = subdf.reset_index()

        coords = {
            'model': list(),
            'fraction': list(),
            'match_rate': list(),
            'phred_mean': list()
        }
        for i in range(0, len(subdf), AUC_STEP):
            sub_df = subdf.loc[i:, :]
            coords['model'].append(self.modelname)
            coords['fraction'].append(len(sub_df)/len(subdf))
            coords['match_rate'].append(np.mean(sub_df['match_rate']))
            coords['phred_mean'].append(np.min(sub_df['phred_mean']))

        coords = pd.DataFrame(coords)
        self.to_csv(coords, reportname)

        return coords

