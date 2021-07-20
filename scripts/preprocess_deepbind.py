import gzip
import csv
#import pandas as pd

seqs = [] # stores dna/rna sequences
descriptors = [] # stores description of sequences
classes = [] # stores 1 for positive and 0 for negative

# parse file of positives
with open(snakemake.input[0], 'r') as fasta_file:
    for line in fasta_file:
        if line[0] == '>':
            descriptors.append(line)
        else:
            seqs.append(line)
            classes.append(1)

# parse file of negatives
with open(snakemake.input[1], 'r') as fasta_file:
    for line in fasta_file:
        if line[0] == '>':
            descriptors.append(line)
        else:
            seqs.append(line)
            classes.append(0)

counter = 0
with open(snakemake.output[0], 'w', newline='') as csvfile:
    # fieldnames = ['FoldID', 'EventID seq', 'Bound']
    # writer = csv.DictWriter(csvfile, fieldnames, delimiter='\t', lineterminator='\n')

    # writer.writeheader()
    # writer.writerow({'FoldID': 'A', 'EventID seq': f'seq_{counter:05d}_peak', 'Bound': (seqs[0] +
    # '\t1')})
    csvfile.write('FoldID\tEventID seq\tBound\t\n')
    for seq, clas in zip(seqs, classes):
        csvfile.write(f'A\tseq_{counter:07d}_peak\t{seq[:100]}\t{clas}\n')
        counter += 1
