import pandas as pd
import gzip

df = pd.read_csv(snakemake.input[0])

# print(df)
# print(snakemake.params)
# print(type(snakemake.params))
# print(snakemake.params[0])


## preprocessing DeepBind
if snakemake.params[0] == 'deepbind':
    print('Writing file for DeepBind...')
    with gzip.open(snakemake.output[0], 'w') as f:
        f.write(b'FoldID\tEventID seq\tBound\t\n')
        for i, seq, clas in zip(df.index, df['seq'], df['class']):
            f.write((f'A\tseq_{i:07d}_peak\t{seq[:101]}\t{clas}\n').encode())
    print('File written!')

## preprocessing iDeepS
elif snakemake.params[0] == 'ideeps':
    print('Writing file for iDeepS...')
    with gzip.open(snakemake.output[0], 'w') as csvfile:
        for descriptor, seq, clas in zip(df['descriptor'], df['seq'], df['class']):
            csvfile.write((descriptor + f'; class:{clas}\n').encode())
            csvfile.write((seq[:101] + '\n').encode())
    print('File written!')

## preprocessing GraphProt2
elif snakemake.params[0] == 'graphprot2':
    print('Writing files for GraphProt2...')
    # needs two files: one with positives and one with negatives
    # no sequences with only N are allowed -> already removed
    # no duplicated chromosome names are allowed -> use index as identifier
    print(snakemake.output['positive'])
    # write positive file
    with open(snakemake.output['positive'], 'w', newline='') as f:
        positives = df[df['class'] == 1]
        for descriptor, seq in zip(positives.index, positives['seq']):
            f.write(f'>{descriptor:05d}\n')
            f.write(seq + '\n')
    print('Positive File written!')

    # write negative file
    with open(snakemake.output['negative'], 'w', newline='') as f:
        negatives = df[df['class'] == 0]
        for descriptor, seq in zip(negatives.index, negatives['seq']):
            f.write('>' + str(descriptor) + '\n')
            f.write(seq + '\n')
    print('Negative file written!')
