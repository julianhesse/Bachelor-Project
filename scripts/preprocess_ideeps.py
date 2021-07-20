import gzip
import csv
#import pandas as pd

seqs = [] # stores dna/rna sequences
descriptors = [] # stores description of sequences
classes = [] # stores 1 for positive and 0 for negative

# parse file of positives
with open('datasets/RBFOX2_HepG2_iDeepS/positive.fasta', 'r') as fasta_file:
    for line in fasta_file:
        if line[0] == '>':
            descriptors.append(line)
        else:
            seqs.append(line)
            classes.append(1)

# parse file of negatives
with open('datasets/RBFOX2_HepG2_iDeepS/negative-1-2.fasta', 'r') as fasta_file:
    for line in fasta_file:
        if line[0] == '>':
            descriptors.append(line)
        else:
            seqs.append(line)
            classes.append(0)

counter = 0
with open('preprocessed/test.seq', 'w', newline='') as csvfile:
    # fieldnames = ['FoldID', 'EventID seq', 'Bound']
    # writer = csv.DictWriter(csvfile, fieldnames, delimiter='\t', lineterminator='\n')

    # writer.writeheader()
    # writer.writerow({'FoldID': 'A', 'EventID seq': f'seq_{counter:05d}_peak', 'Bound': (seqs[0] +
    # '\t1')})
    for descriptor, seq, clas in zip(descriptors, seqs, classes):
        csvfile.write(descriptor[:-1] + f'; class:{clas}\n')
        csvfile.write(seq[:101] + '\n')
