import matplotlib.pyplot as plt
import pandas as pd

cell_line_0 = snakemake.input['cell_line_0']
cell_line_1 = snakemake.input['cell_line_1']
method = snakemake.wildcards['method']
mode = snakemake.wildcards['mode']


## prepare data ##

df_0 = pd.read_csv(cell_line_0).set_index(['RBP', 'fold'])
df_0 = df_0.loc[(slice(None), 'mean'), :]

df_1 = pd.read_csv(cell_line_1).set_index(['RBP', 'fold'])
df_1 = df_1.loc[(slice(None), 'mean'), :]

result = pd.merge (df_0, df_1, on=['RBP'])

plt.plot([0,1], [0,1], c="orange")
# when merging methods of df_0 get '_x' appended
# when merging methods of df_1 get '_y' appended
plt.scatter(result[method + '_x'], result[method + '_y'])
plt.xlim(0.5, 1.)
plt.ylim(0.5, 1.)

if mode == 'roc_auc':
    plt.xlabel(f"Mean auROC over 5 folds: {method}, {snakemake.wildcards.cell_line_0}")
    plt.ylabel(f"Mean auROC over 5 folds: {method}, {snakemake.wildcards.cell_line_1}")
else:
    plt.xlabel(f"Mean AP over 5 folds: {method}, {snakemake.wildcards.cell_line_0}")
    plt.ylabel(f"Mean AP over 5 folds: {method}, {snakemake.wildcards.cell_line_1}")

plt.savefig(snakemake.output[0], dpi=300)
plt.show()
