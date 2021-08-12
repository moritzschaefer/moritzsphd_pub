import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


class SharedColorPlot:
    def __init__(self, ncols, nrows, num_colors=20):
        self.fig, self.axes = plt.subplots(ncols, nrows, squeeze=False)
        self.cluster_color_order = []
        self.cmap = plt.get_cmap('jet', num_colors)

    def plot(self,
             y,
             x,
             df=None,
             kind='pie',
             ylabel='',
             labelpad=None,
             proportion=None,
             title='',
             num_labels=4):
        if df is None:
            self.fig.delaxes(self.axes[y, x])
            return

        colors = []
        for cluster, rest in df.iterrows():
            if cluster not in self.cluster_color_order:
                self.cluster_color_order.append(cluster)

            colors.append(self.cmap(self.cluster_color_order.index(cluster)))

        # TODO I don't know if this works
        if proportion and proportion < 1:
            filling = (df['count'].sum() / proportion) - df['count'].sum()
            df.append(pd.Series({'count': filling}, name=''))

        ax = df.plot(
            kind=kind,
            ax=self.axes[y][x],
            y=0,
            legend=False,
            rotatelabels=False,
            labels=df.index[:num_labels].tolist() +
            [''] * (len(df) - (num_labels + 1)) + df.index[-1:].tolist(),
            colors=colors)

        ax.set_ylabel(ylabel, labelpad=labelpad)
        ax.set_title(title)

    def legend(self, **kwargs):
        patches = [
            mpatches.Patch(
                color=self.cmap(self.cluster_color_order.index(cluster)),
                label=cluster) for cluster in self.cluster_color_order
        ]
        self.fig.legend(
            patches,
            self.cluster_color_order,
            loc='upper right',
            frameon=True,
            bbox_to_anchor=(0.95, 0.95),
            borderaxespad=0,
            **kwargs)
