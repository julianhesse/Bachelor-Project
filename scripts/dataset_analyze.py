import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

rbp_files = dict(snakemake.input)

values = []

for rbp, f_name in rbp_files.items():
    df = pd.read_csv(f_name)

    positives = len(df[df['class'] == 1])
    negatives = len(df[df['class'] == 0])

    neg_ratio = negatives / (positives + negatives)

    values.append((positives, negatives, neg_ratio))


values = np.array(values)
plt.plot([0,np.max(values)], [0,np.max(values)], c="orange", zorder=1)
plt.scatter(values[:, 0], values[:, 1], s=8, alpha=0.5, zorder=2)
plt.xlabel("number of positives")
plt.ylabel("number of negatives")
plt.title(f'{snakemake.wildcards["cell_line"]}')

plt.savefig(snakemake.output[0], dpi=300)
