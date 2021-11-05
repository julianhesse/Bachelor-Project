import pandas as pd
import numpy as np

rbp_files = dict(snakemake.input)

values = {}

for rbp, f_name in rbp_files.items():
    df = pd.read_csv(f_name)

    positives = len(df[df['class'] == 1])
    negatives = len(df[df['class'] == 0])

    neg_ratio = negatives / (positives + negatives)

    values[rbp] = [positives, negatives, neg_ratio]


df = pd.DataFrame(values, index=['num_pos', 'num_neg', 'neg_ratio']).T

df.to_csv(snakemake.output[0])
