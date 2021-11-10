import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

methods = snakemake.config['methods']
boundary = int(snakemake.wildcards['boundary'])/10

df = pd.read_csv(snakemake.input[0])

for method in methods:
    df.loc[df[method] < boundary, method] = 0
    df.loc[df[method] > 0, method] = 1

df = df.astype('int')


arr = np.array([
    df['class'],
    *[df[method] for method in methods]
])

fig = plt.figure(figsize=(8,6))
plt.matshow(arr, aspect="auto")
plt.yticks(np.arange(0, len(methods)+1),['label', *methods])
fig.suptitle(f'')
plt.savefig(snakemake.output[0], dpi=300)
