import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import transforms
from scikitplot.decomposition import plot_pca_2d_projection
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def maplot(data, mean='baseMean', diff='log2FoldChange', p='padj', p_threshold=0.05, ax=None, palette={'up': 'red', 'down': 'blue', 'unchanged': 'black'}, label_ss=False):
    data = data.copy()
    data['hue'] = 'unchanged'
    data.loc[(data[p] < p_threshold) & (data[diff] > 0.5), 'hue'] = 'up'
    data.loc[(data[p] < p_threshold) & (data[diff] < -0.5), 'hue'] = 'down'
    data['log2(mean expression + 1)'] = np.log2(data[mean] + 1)
    volcanoplot(data=data, ax=ax, show_n='y', x='log2(mean expression + 1)', y=diff, s=7, palette=palette, label_ss=label_ss)


def volcanoplot(data, x='log2FoldChange', y='-log10(adjusted pvalue)', show_n=False, ax=None, alpha=1.0, hue=None, palette={'up': 'red', 'down': 'blue', 'unchanged': 'black'}, y_threshold=-np.log10(0.05), x_threshold=0.5, s=40, label_ss=False):
    '''
    Given a df, plot a volcano plot.

    Either provide a field for the coloring on your own (hue parameter), or provide a threshold via y_threshold and x_threshold to compute the color

    Default fields: 'log2FoldChange' and '-log10(adjusted pvalue)'
    cutoffs for p-value is 0.05 and for log2FoldChange 0.5
    '''
    data = data.copy()
    if not hue:
        hue = 'hue'
    if hue not in data.columns:
        data[hue] = 'unchanged'
        data.loc[(data[x] > x_threshold) & (data[y] > y_threshold), hue] = 'up'
        data.loc[(data[x] < -x_threshold) & (data[y] > y_threshold), hue] = 'down'

    # sort by hue so that insignificant is drawn first
    data['hue'] = pd.Categorical(data['hue'], categories=['unchanged', 'up', 'down'], ordered=True)
    data = data.sort_values('hue')

    ax = sns.scatterplot(data=data, x=x, y=y, hue=hue, palette=palette, s=s, edgecolor="none",
                         legend=False, alpha=alpha, ax=ax, rasterized=True, linewidth=0)
    if show_n:
        if show_n == 'y':
            x_up, x_down = data[x].max(), data[x].max()
            y_up = min(data[y].max() * 0.8, ax.get_ylim()[1])
            y_down = max(data[y].min() * 0.8, ax.get_ylim()[0])
            ha_up = 'right'
            ha_down = 'right'
        else:
            y_up, y_down = data[y].max(), data[y].max()
            x_up = min(data[x].max() * 0.9, ax.get_xlim()[1])
            x_down = max(data[x].min() * 0.9 , ax.get_xlim()[0])
            ha_up = 'right'
            ha_down = 'left'

        ax.text(x_down, y_down, f'n={(data[hue] == "down").sum()}', color=palette['down'], horizontalalignment=ha_down)
        ax.text(x_up, y_up, f'n={(data[hue] == "up").sum()}', color=palette['up'], horizontalalignment=ha_up)

    if label_ss:
        for index, row in data[~(data[hue] == 'unchanged')].iterrows():
            ax.text(row[x], row[y], index)

    sns.despine()

    return ax


def pcaplot(df, dimensions=[0, 1], scaling=None, ax=None, label_dots=True, cmap='Spectral', stripext=False):
    '''
    Requires a normalized library in the first place (i.e. each sample should be CPM (or so) normalized)
    # Normalization for PCA: https://chanzuckerberg.github.io/scRNA-python-workshop/preprocessing/02-normalization.html

    scaling: - 'standard': use standard scaler to scale the features
             - 'log2': use log2 eof all values
             - 'log2-standard' First apply log2, then standardize
             - None: don't scale
             'standard' and 'log2' help to avoid PCA domination by large genes

    '''
    if len(dimensions) != 2:
        raise ValueError('Only 2D plot supported')

    pca_data = df.T.values

    if scaling:
        if 'log2' in scaling:
            pca_data = np.log2(pca_data + 1 )

        if 'standard' in scaling:
            scaler = StandardScaler()
            pca_data = scaler.fit_transform(pca_data)

    # df.index.name = 'miRNA'
    pca = PCA()
    pca.fit(pca_data)

    explained_variance = np.sum(pca.explained_variance_ratio_[dimensions])

    if len(dimensions) == 2:  #
        if stripext:
            labels = np.array([s[:s.rfind('_')] for s in df.columns])
        else:
            labels = df.columns
        ax = plot_pca_2d_projection(pca, pca_data, labels, ax=ax, label_dots=label_dots, cmap=cmap, feature_labels=None, dimensions=dimensions)
    else:
        ax = plot_pca_3d_projection(pca, pca_data, df.columns, ax=ax, label_dots=label_dots, cmap=cmap, feature_labels=None)

    ax.set_title(f'Explained variance: {explained_variance}')
    return ax, pca.explained_variance_ratio_[dimensions]


def grouped_barplot(data, x, hue, y, yerr, split_yaxis=None, ax=None, colors=None, group_offsets=None, group_labels=None):
    '''
    This is a desaster! This function was not designed for what it is used now. Now it's used for qpcranalysis plots..

    Bar plot function with custom error bars
    :split_yaxis: split y_axis given a key/column to group by
    :group_offsets: make a spacer between bars at given positions
    '''
    if split_yaxis and ax:
        raise ValueError('Cannot use ax and split_yaxis at the same time')

    def _single_grouped_barplot(data, x, hue, y, yerr, ax):
        u = data[x].unique()
        x = np.arange(len(u))
        subx = data[hue].unique()
        offsets = (np.arange(len(subx))-np.arange(len(subx)).mean())/(len(subx)+1.)
        width = np.diff(offsets).mean()
        if group_offsets is not None:
            for goffset in group_offsets:
                offsets += np.array([width * 0.5 if i >= goffset else 0.0 for i in range(len(offsets))])
        for i, gr in enumerate(subx):
            dfg = data[data[hue] == gr]
            ax.bar(x+offsets[i], dfg[y].values, width=width,
                   label="{} {}".format(hue, gr), yerr=dfg[yerr].values,
                   color=colors[i] if colors is not None else None)
        # ax.set_xlabel(x)
        # ax.set_ylabel(y)
        if len(u) == 1:
            ax.set_title(u[0])
            ax.set_xticks(x+offsets)
            ax.set_xticklabels(subx, rotation=50, ha='right')
        else:
            ax.set_xticks(x)
            ax.set_xticklabels(u)
        if group_labels is not None:
            ax.set_xticks(offsets[[go - 1 for go in group_offsets]])
            ax.set_xticklabels(group_labels, rotation=20, ha='right')
        ax.set_ylim([0, data[y].max() * 1.1])

    if split_yaxis and len(data.groupby(split_yaxis)) > 1:
        groups = data.groupby(split_yaxis)
        fig, axes = plt.subplots(1, len(groups), figsize=(3.4 * len(groups), 4))
        for ax, group in zip(axes, groups):
            _single_grouped_barplot(group[1].reset_index(), x, hue, y, yerr, ax)
    else:
        if not ax:
            fig, ax = plt.subplots()
        _single_grouped_barplot(data, x, hue, y, yerr, ax)
    return fig


def significance_level(ax, x1, x2, level, increase_plot_y=True):
    '''
    Add significance level to pyplot.Axis
    The y_lim of the axis is changed such that subsequent calls to the function automatically increase in height

    :ax: The pyplot.Axis to draw in
    :x1: Start x position of the bar
    :x2: End x position of the bar
    :level: significance level
    '''
    yticks = ax.get_yticks()

    y_min, y_max = ax.get_ylim()
    fig_width = ax.get_figure().get_figwidth()
    fig_height = ax.get_figure().get_figheight()

    ax.annotate("", xy=(x1, y_max), xycoords='data',
                xytext=(x2, y_max), textcoords='data',
                arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                                connectionstyle=f"bar,fraction={1 / ((x2-x1) * fig_width) }"))
    text_y = y_max + (y_max - y_min) * 0.5 / fig_height
    ax.text(x1 + (x2-x1)/2, text_y, f'p < {level:.4f}',
            horizontalalignment='center',
            verticalalignment='center')

    if increase_plot_y:
        ax.set_ylim([y_min, y_max + (y_max - y_min) * 0.13])

    ax.set_yticks(yticks)

def colored_text(x, y, strings, colors, ax=None, **kwargs):
    """
    Take a list of *strings* and *colors* and place them next to each
    other, with text strings[i] being shown in colors[i].

    Parameters
    ----------
    x, y : float
        Text position in data coordinates.
    strings : list of str
        The strings to draw.
    colors : list of color
        The colors to use.
    ax : Axes, optional
        The Axes to draw into. If None, the current axes will be used.
    **kwargs
        All other keyword arguments are passed to plt.text(), so you can
        set the font size, family, etc.
    """
    if ax is None:
        ax = plt.gca()
    t = ax.transData
    canvas = ax.figure.canvas

    for s, c in zip(strings[::-1], colors[::-1]):
        text = ax.text(x, y, s + " ", color=c, transform=t, **kwargs)

        # Need to draw to update the text position.
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()
        t = transforms.offset_copy(
            text.get_transform(), y=ex.height/2, units='dots')
