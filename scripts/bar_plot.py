import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

file = snakemake.input[0]

folds = snakemake.config['folds']
cell_line = snakemake.wildcards['cell_line']
methods = snakemake.config['methods']
mode = snakemake.wildcards['mode']

## prepare data ##

df = pd.read_csv(file).set_index(['RBP', 'fold'])

df = df.loc[(slice(None), 'mean'), :]

df = df.sort_values(by=methods[0], ascending=False)

## plot ##

labels = list(df.index.get_level_values(0))

x = np.arange(len(labels)) * 1.5
width = 0.4

fig, ax = plt.subplots()

bars = []
for i, method in enumerate(methods):
    x_moved = x + (i) * width - (len(methods)-1)/2 * width
    bar = ax.bar(x_moved, df[method], width, label=method)
    bars.append(bar)


if mode == 'roc_auc':
    ax.set_ylabel('mean AUC ROC for each dataset')
else:
    ax.set_ylabel('mean AP for each dataset')
ax.set_xlabel('methods')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
ax.set_axisbelow(True)

for tick in ax.get_xticklabels():
    tick.set_rotation(70)

start = 0.4
ax.set_ylim([start,1.])
ax.set_yticks(np.arange(start,1.1, step=0.1))
# ax.yaxis.grid(True)

# for bar in bars:
#     ax.bar_label(bar, padding=3, fmt='%.2f')

fig.tight_layout()
plt.grid(axis='y', zorder=0)
plt.show()
fig.savefig(snakemake.output[0], dpi=400)
