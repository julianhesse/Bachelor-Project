import pandas as pd
import numpy as np
import logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
        format='%(asctime)s %(name)s - %(levelname)s: %(message)s')

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
logging.info('Positives parced...')

# parse file of negatives
with open(snakemake.input[1], 'r') as fasta_file:
    for line in fasta_file:
        if line[0] == '>':
            descriptors.append(line[:-1])
        else:
            seqs.append(line[:-1])
            classes.append(0)
logging.info('Negatives parced...')

df = pd.DataFrame(data={
    'seq': seqs,
    'descriptor': descriptors,
    'class': classes
})

# remove sequences with only N
df = df[df['seq'].str.match(r'^N+N$') == False]
if sum(df['seq'].str.match(r'N+')) > 0:
    logging.info('Sequences with N remaining:')
    logging.info(df['seq'].str.match(r'^N+N$'))
else:
    logging.info('Sequences with only N removed...')

# remove duplicates
num_dup = sum(df.duplicated())
num_dup_chrom = sum(df.duplicated(subset=['descriptor']))

logging.info('Number of duplicates: %d', num_dup)
logging.info('Number of duplicates of chromosomes: %d', num_dup_chrom)

if snakemake.config['drop_duplicates']:
    # only drop duplicate sequences with the same class
    logging.info('Dropped duplicates...')
    df = df.drop_duplicates()
else:
    logging.info('No duplicates dropped...')

if snakemake.config['drop_duplicate_chromosomes']:
    # drop all duplicate sequences
    # keeps positive versions of the chromosome
    # because positives are first in the dataframe
    logging.info('Dropped duplicate chromosomes')
    df = df.drop_duplicates(subset=['descriptor'])
else:
    logging.info('No duplicate chromosomes dropped...')

df.reset_index(drop=True, inplace=True)

# create folds
# perserve ratio of classes in folds

folds = snakemake.params['folds']
logging.info(f'Partition dataset into {folds} folds ...')

ids_pos = df[df["class"] == 1].index.to_numpy()
ids_neg = df[df["class"] == 0].index.to_numpy()

np.random.seed(snakemake.params['seed'])
ids_pos = np.random.permutation(ids_pos)
ids_neg = np.random.permutation(ids_neg)

folds_pos = np.array_split(ids_pos, folds)
folds_neg = np.array_split(ids_neg, folds)

folds_ids = [np.concatenate((pos,neg)) for pos, neg in zip(folds_pos, folds_neg)]

fold_column = np.zeros(len(df), dtype=int)

for i, ids in enumerate(folds_ids):
    fold_column[ids] = i

logging.info("Size of folds:")
unique, counts = np.unique(fold_column, return_counts=True)
logging.info(dict(zip(unique, counts)))

df['fold'] = fold_column

logging.info(df)

logging.info('Folds are defined...')


# write remaining entries to csv
df.to_csv(snakemake.output[0])
logging.info('CSV file created!')
logging.info('Dataset preprocessed!')
