import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyperclip
import seaborn as sns


class QPCRAnalysis:
    def __init__(self,
                 data,
                 genes,
                 samples,
                 ntc_cols=True,
                 columns_per_sample=2):
        '''
        data: can be a filename to a tsv as exported by the lightcycler or a dataframe
        '''
        self.samples = samples
        self.genes = genes
        if isinstance(data, pd.DataFrame):
            df = data
        else:
            with open(data, 'r') as f:
                self.experiment_name = re.match('Experiment: ([^ ]*)',
                                                f.readline()).groups()[0]
            df = pd.read_csv(data, sep='\t', skiprows=1)

        df.drop('Name', inplace=True, axis=1)

        df['ypos'] = df['Pos'].apply(lambda pos: ord(pos[0]) - ord('A'))
        df['ypos'] -= df['ypos'].min()

        ntc_row = df['ypos'].max()

        df['xpos'] = df['Pos'].apply(lambda pos: int(pos[1:]))
        # only consider non-NTCs here
        df['xpos'] -= df.loc[df['ypos'] < ntc_row, 'xpos'].min()

        if ntc_cols:
            if ntc_cols is True:
                ntc_cols = list(range(0, len(genes) * 2, 2))

            self.ntcs = pd.Series(
                df.loc[df['xpos'].isin(ntc_cols) & (df['ypos'] == ntc_row),
                       'Cp'].tolist(),
                index=genes)
        else:
            self.ntcs = None
        df = df[(df.ypos < ntc_row)
                & (df.xpos < len(samples) * columns_per_sample)]

        df['gene'] = self._assign_genes(df, genes, columns_per_sample)

        df['sample'] = pd.Categorical(
            df.xpos.apply(lambda xpos: samples[xpos // columns_per_sample]),
            categories=samples,
            ordered=True)

        self.df = df[['sample', 'gene', 'Cp', 'xpos', 'ypos']]

    def dropna(self):
        na_index = self.df.index[self.df.Cp.isna()]
        print(f'Deleting {len(na_index)} Nans')
        self.df.drop(na_index, inplace=True)

    def _assign_genes(self, df, genes, columns_per_sample):
        def _gene(row):
            index = ((row['ypos'] // 3) * columns_per_sample +
                     (row['xpos'] % columns_per_sample))
            try:
                return genes[index]
            except IndexError:
                print(index)
                raise

        return df.apply(_gene, axis=1)

    def outliers(self):
        df = self.df.copy()
        triplet_deviations = df.groupby(['sample', 'gene'])['Cp'].std()

        df['index'] = df.index
        outlier_triplets = df.set_index(
            ['sample',
             'gene']).loc[triplet_deviations.index[triplet_deviations > 0.5]]
        outlier_ids = []
        for name, group in outlier_triplets.groupby(level=['sample', 'gene']):
            diffs = (group.Cp - group.Cp.mean())
            if len(diffs) == 3:
                sorted_diffs = diffs.abs().reset_index(drop=True).sort_values(ascending=False)
                if sorted_diffs.iloc[0] > 1.5 * sorted_diffs.iloc[1]:
                    index = group.iloc[sorted_diffs.index[0]]['index']
                    outlier_ids.append(int(index))


                # print(diffs.iloc[lid], df[df['index'] == index])
            # elif len(diffs) == 2:  # does this make sense? if we have two, then leave the two!
            #     # delete the larger one
            #     index = group.iloc[group.Cp.reset_index(
            #         drop=True).argmax()]['index']

        return outlier_ids

    def plot_outliers(self):  # TODO
        triplet_deviations = self.df.groupby(['sample', 'gene']).std()['Cp']

        outlier_triplets = self.df.set_index(['sample', 'gene']).loc[
            triplet_deviations.index[triplet_deviations > 0.5]].reset_index()
        outlier_triplets['triplet_name'] = outlier_triplets['sample'].str.cat(
            outlier_triplets.gene, sep='_')
        sns.swarmplot(data=outlier_triplets.reset_index(),
                      x='triplet_name',
                      y='Cp')
        plt.xticks(rotation=90)

    def plot_cps(self, exclude_genes=[]):
        plt.subplots(figsize=(10, 7))
        # g = sns.barplot(y='Cp', hue='sample', x='gene', data=self.df.sort_values('sample'), dodge=True)
        plot_df = self.df.sort_values('sample')
        plot_df = plot_df[~plot_df.gene.isin(exclude_genes)]
        g = sns.barplot(y='Cp',
                        hue='sample',
                        x='gene',
                        data=plot_df,
                        dodge=True)
        plt.xticks(rotation=30)
        plt.title('Raw Cp values')
        plt.legend(ncol=3, loc='upper center', bbox_to_anchor=[0.5, -0.15])

    def normalized_df(self, normalizers, exclude_genes=[], exclude_samples=[], include_samples=None, include_genes=None, norm_to_one=None):
        '''
        :include_samples: takes precedence over exclude_samples. Can be used to determ̀ine the plotting order
        :include_genes: takes precedence over exclude_genes. Can be used to determ̀ine the plotting order
        '''
        assert len(normalizers) > 0, 'At least one normalizer is required'
        means = self.df.groupby(['sample', 'gene']).mean().unstack()['Cp']
        stds = self.df.groupby(['sample', 'gene']).std().unstack()['Cp']

        summed_normalizers = None
        for normalizer in normalizers:
            if summed_normalizers is None:
                summed_normalizers = means[normalizer]
            else:
                summed_normalizers = summed_normalizers + means[normalizer]
        summed_normalizers = summed_normalizers / len(normalizers)

        delta = means.subtract(
            summed_normalizers, axis=0
        )  # normalization like this results in geometric mean normalization on mRNA expression level
        non_normalization_genes = [
            gene for gene in self.genes if gene not in normalizers
        ]
        if include_genes:
            plot_genes = [
                gene for gene in include_genes
                if gene in non_normalization_genes
            ]
        else:
            plot_genes = [
                gene for gene in non_normalization_genes
                if gene not in exclude_genes
            ]
        delta_low = delta - stds
        delta_high = delta + stds
        q = np.power(2, -delta)
        q_low = np.power(2, -delta_low)
        q_high = np.power(2, -delta_high)
        q_std = pd.concat(dict(q=q, q_low=q_low, q_high=q_high),
                          axis=1).std(axis=1, level=1)
        # fig, ax = plt.subplots(figsize=(10, 5))
        plot_df = q[plot_genes].stack().reset_index().rename(
            columns={0: 'expression'})
        plot_df['error'] = q_std[plot_genes].stack().reset_index()[0]

        if include_samples is not None:
            plot_df = plot_df[plot_df['sample'].isin(include_samples)]
            plot_df['sample'] = pd.Categorical(plot_df['sample'],
                                               categories=include_samples,
                                               ordered=True)
            plot_df.sort_values('sample', inplace=True)
        else:
            plot_df = plot_df[~plot_df['sample'].isin(exclude_samples)]

        if include_genes:
            plot_df['gene'] = pd.Categorical(plot_df['gene'],
                                               categories=plot_genes,
                                               ordered=True)
            plot_df.sort_values('gene', inplace=True, kind='mergesort')  # mergesort is stable


        if norm_to_one:
            normed = plot_df.set_index(['sample', 'gene']).groupby('gene').apply(lambda group: pd.DataFrame({
                'expression': group['expression'] / group.loc[norm_to_one, 'expression'],
                'error': group['error'] / group.loc[norm_to_one, 'expression']}))
            return normed.reset_index()
        else:
            return plot_df

    def plot_normalized(self,
                        normalizers,
                        exclude_genes=[],
                        exclude_samples=[],
                        include_samples=None,
                        include_genes=None,
                        colors=None,
                        legend=True,
                        norm_to_one=None,
                        **kwargs
    ):
        '''
        This has been built using the Ciaudo Lab Excel qPCR analysis sheet as template
        '''
        from moritzsphd.plot import grouped_barplot
        plot_df = self.normalized_df(normalizers, exclude_genes, exclude_samples, include_samples, include_genes, norm_to_one)

        fig = grouped_barplot(data=plot_df,
                              y='expression',
                              x='gene',
                              hue='sample',
                              yerr='error',
                              split_yaxis='gene',
                              colors=colors,
                              **kwargs
        )
        #sns.barplot(data=plot_df, y='expression', x='gene', hue='sample')
        #plt.errorbar(x=plot_df['gene'], y=plot_df['expression'], fmt='none', yerror=plot_df['error'], ecolor='k', elinewidth=2)
        plt.subplots_adjust(top=0.85)
        plt.suptitle(
            f'gene expression normalized by {"geometric mean of " if len(normalizers) > 1 else ""}{" and ".join(normalizers)}'
        )
        sns.despine()
        if legend:
            plt.legend(ncol=3, loc='upper center', bbox_to_anchor=[0.5, -0.1])
        return fig


    def drop_outliers(self):
        self.df.drop(self.outliers(), inplace=True)

    def plot_heatmap(self):
        fig, ax = plt.subplots(figsize=(15, 7))
        plot_df = self.df[['Cp', 'xpos', 'ypos']].pivot(index='ypos',
                                                        columns='xpos')
        sns.heatmap(plot_df)

    def export_excel(self, normalizer_gene, readout_gene):
        df = self.df.copy()
        normalizer = df.loc[df.gene == normalizer_gene].set_index(
            'sample').sort_index().Cp
        xx = df.loc[df.gene == readout_gene].set_index('sample').sort_index()
        values = []
        for i, x in enumerate(xx.iterrows()):
            values.append(x[1].Cp)
            if i % 3 == 2:  # now gapdh
                for j in range(i - 2, i + 1):
                    values.append(normalizer.iloc[j])

        if len(values) == 0:
            raise ValueError(f'gene {readout_gene} does not exist')
        pyperclip.copy('\n'.join(
            ['' if pd.isnull(v) else str(v) for v in values]))

    def normalized_prism_df(self, grouping_func, normalizer_genes=None, control_group='WT', repl_normalization=True):
        '''
        '''
        df = self.normalized_df(normalizer_genes)
        df['grouping_series'] = df.apply(grouping_func, axis=1)

        # bring into PRISM form
        df = df.groupby(['grouping_series', 'gene'])['expression'].apply(lambda v: pd.Series(list(v))).unstack(-2).unstack(-1)

        # normalize to WT
        if repl_normalization:
            df = df / df.loc[control_group]
        else:
            wt_mean = df.loc[control_group].unstack(-1).mean(axis=1)
            df = df.groupby(axis=1, level=0).apply(lambda v: v/wt_mean[v.name])

        return df

    def raw_prism_df(self):
        df = self.df.copy()
        df['repl'] = df['ypos'] % 3

        return df.set_index(['sample', 'gene', 'repl'])['Cp'].unstack(-2).unstack(-1)
