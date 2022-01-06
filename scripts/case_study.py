import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn2

methods = snakemake.config['methods']


# load data
df = pd.read_csv(snakemake.input[0])
folds = snakemake.config['folds']

for method in methods:
    df.loc[df[method] < 0.5, method] = 0
    df.loc[df[method] > 0, method] = 1

df = df.astype('int')


## all methods correct ##


def get_sets(df):
    classes = df['class']
    sets = []

    for method in methods:
        sets.append(set(df[df[method] == classes].index))

    return sets

fig, axs = plt.subplots(folds+1, 2, figsize=(4*2, 4*(folds)))

df_tmp = df[df['class'] == 0]

sets = get_sets(df_tmp)

if len(methods) == 3:
    venn3(sets, methods, ax=axs[0,0])
else:
    venn2(sets, methods, ax=axs[0,0])

axs[0,0].set_title('negatives')

df_tmp = df[df['class'] == 1]

sets = get_sets(df_tmp)

if len(methods) == 3:
    venn3(sets, methods, ax=axs[0,1])
else:
    venn2(sets, methods, ax=axs[0,1])

axs[0,1].set_title('positives')

for i in range(folds):
    df_tmp = df[(df['class'] == 0) & (df['fold'] == i)]

    sets = get_sets(df_tmp)

    if len(methods) == 3:
        venn3(sets, methods, ax=axs[i+1,0])
    else:
        venn2(sets, methods, ax=axs[i+1,0])

    df_tmp = df[(df['class'] == 1) & (df['fold'] == i)]

    sets = get_sets(df_tmp)

    if len(methods) == 3:
        venn3(sets, methods, ax=axs[i+1,1])
    else:
        venn2(sets, methods, ax=axs[i+1,1])

fig.suptitle(f'{snakemake.wildcards.dataset} - venn-diagrams for correct predictions')
fig.savefig(snakemake.output[0], dpi=300)
plt.show()


if len(methods) == 3:
    venn3(get_sets(df), methods)
else:
    venn2(get_sets(df), methods)
plt.savefig(snakemake.output[1], dpi=300)
plt.show()
