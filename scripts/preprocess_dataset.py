import pandas as pd
import numpy as np

seqs = [] # stores dna/rna sequences
descriptors = [] # stores description of sequences
classes = [] # stores 1 for positive and 0 for negative

# parse file of positives
with open(snakemake.input[0], 'r') as fasta_file:
    for line in fasta_file:
        if line[0] == '>':
            descriptors.append(line[:-1])
        else:
            seqs.append(line[:-1])
            classes.append(1)
print('Positives parced...')

# parse file of negatives
with open(snakemake.input[1], 'r') as fasta_file:
    for line in fasta_file:
        if line[0] == '>':
            descriptors.append(line[:-1])
        else:
            seqs.append(line[:-1])
            classes.append(0)
print('Negatives parced...')

df = pd.DataFrame(data={
    'seq': seqs,
    'descriptor': descriptors,
    'class': classes
})

# remove sequences with only N
df = df[df['seq'].str.match(r'^N+N$') == False]
if sum(df['seq'].str.match(r'N+')) > 0:
    print('Sequences with N remaining:')
    print(df['seq'].str.match(r'^N+N$'))
    print()
else:
    print('Sequences with only N removed...')

# remove duplicates
num_dup = sum(df.duplicated())
num_dup_chrom = sum(df.duplicated(subset=['descriptor']))

print('Number of duplicates:', num_dup)
print('Number of duplicates of chromosomes:', num_dup_chrom)

if snakemake.config['drop_duplicates']:
    # only drop duplicate sequences with the same class
    print('Dropped duplicates...')
    df = df.drop_duplicates()
else:
    print('No duplicates dropped...')

if snakemake.config['drop_duplicate_chromosomes']:
    # drop all duplicate sequences
    # keeps positive versions of the chromosome
    # because positives are first in the dataframe
    print('Dropped duplicate chromosomes')
    df = df.drop_duplicates(subset=['descriptor'])
else:
    print('No duplicate chromosomes dropped...')

df.reset_index(inplace=True)

# create split: test and train

print('Splitting dataset ...')
split = df.groupby(['class']).sample(frac=snakemake.params['frac'], random_state=snakemake.params['seed'])
part_class_0 = sum(split['class']==0)/sum(df['class'] == 0)
part_class_1 = sum(split['class']==1)/sum(df['class'] == 1)
# make sure negatives and positives are split independently
assert(abs(snakemake.params['frac'] - part_class_0) < 0.01, "Negatives and positives not split independently!")
assert(abs(snakemake.params['frac'] - part_class_1) < 0.01, "Negatives and positives not split independently!")

msk_split = np.zeros(len(df.index), dtype=bool)
msk_split[split.index] = True
df['test'] = msk_split

print(df)

print('Train and test dataset defined')


# write remaining entries to csv
df.to_csv(snakemake.output[0])
print('CSV file created!')
