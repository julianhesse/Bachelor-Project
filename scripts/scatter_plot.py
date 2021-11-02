import matplotlib.pyplot as plt
import pandas as pd

file = snakemake.input[0]

folds = snakemake.config['folds']
cell_line = snakemake.wildcards['cell_line']
method_0 = snakemake.wildcards['method_0']
method_1 = snakemake.wildcards['method_1']
mode = snakemake.wildcards['mode']


## prepare data ##

df = pd.read_csv(file).set_index(['RBP', 'fold'])

df = df.loc[(slice(None), 'mean'), :]

plt.plot([0,1], [0,1], c="orange")
plt.scatter(df[method_0], df[method_1])
plt.xlim(0.5, 1.)
plt.ylim(0.5, 1.)

if mode == 'roc_auc':
    plt.xlabel(f"Mean auROC over 5 folds: {method_0}")
    plt.ylabel(f"Mean auROC over 5 folds: {method_1}")
else:
    plt.xlabel(f"Mean AP over 5 folds: {method_0}")
    plt.ylabel(f"Mean AP over 5 folds: {method_1}")

plt.savefig(snakemake.output[0], dpi=300)
plt.show()
