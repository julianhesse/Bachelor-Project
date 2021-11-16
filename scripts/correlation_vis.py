import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

dataset_stats = snakemake.input['stats']
dataset_performance = snakemake.input['performance']
statistic = snakemake.wildcards['statistic']
mode = snakemake.wildcards['mode']
cell_line = snakemake.wildcards['cell_line']
methods = snakemake.config['methods']

df_stats = pd.read_csv(dataset_stats, index_col=0).sort_index()
df_performance = pd.read_csv(dataset_performance).set_index(['RBP', 'fold'])
df_performance = df_performance.loc[(slice(None), statistic), :]

size = df_stats['num_pos'] + df_stats['num_neg']

scale=snakemake.params['scale']
fig, axs = plt.subplots(1,len(methods), figsize=(scale*len(methods),scale*1), sharey='row')

# axs = list(axs)

for i, method in enumerate(methods):
    # axs[i].plot([0,1], [0,1], c="orange", zorder=1)


    ## scatter with points for each dataset ##
    ax = axs[i]
    ax.scatter(size, df_performance[method], s=8, alpha=0.5, label='dataset')

    ## calculate a regression line lstsq ##
    A = np.vstack([size, np.ones(len(size))]).T
    m, b = np.linalg.lstsq(A, df_performance[method], rcond=None)[0]
    ax.plot(size, m*size + b, c='orange',  label='fitted line')

    # p = np.polynomial.polynomial.Polynomial.fit(size, df_performance[method], 2)
    # x = np.linspace(0, np.max(size), 30)
    # ax.plot(x, p(x), c='red',  label='fitted line 2')

    ax.set_title(f'{method}')
    ax.set_xlabel(f'Size of dataset from {cell_line}')
    ax.set_ylabel(f'{statistic} {mode} {method}')
    ax.legend()


fig.tight_layout()
plt.savefig(snakemake.output[0], dpi=300)
