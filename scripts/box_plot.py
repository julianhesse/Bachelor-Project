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

## plot ##

fig, ax = plt.subplots()

boxplot = ax.boxplot(df,
           vert=True,
           patch_artist=True,
           labels=list(df.columns))

ax.set_ylim([0,1.])
ax.set_yticks(np.arange(0,1.1, step=0.1))
if mode == 'roc_auc':
    plt.xlabel(f"Mean auROC over 5 folds")
    plt.ylabel(f"Mean auROC over 5 folds")
else:
    plt.xlabel(f"Mean AP over 5 folds")
    plt.ylabel(f"Mean AP over 5 folds")
ax.set_xlabel('methods')
ax.yaxis.grid(True)

# fill with colors
colors = ['pink', 'lightblue', 'lightgreen']
for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)

fig.savefig(snakemake.output[0], dpi=300)
