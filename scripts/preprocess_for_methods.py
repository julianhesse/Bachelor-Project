import pandas as pd
import gzip

df = pd.read_csv(snakemake.input[0])

# print(df)
# print(snakemake.params)
# print(type(snakemake.params))
# print(snakemake.params[0])

## split into train and test data


## preprocessing DeepBind
if snakemake.params[0] == 'deepbind':
    print('Writing files for DeepBind...')
    data = df[df['test'] == False]
    with gzip.open(snakemake.output['train'], 'w') as f:
        f.write(b'FoldID\tEventID seq\tBound\t\n')
        for i, seq, clas in zip(data.index, data['seq'], data['class']):
            f.write((f'A\tseq_{i:07d}_peak\t{seq[:101]}\t{clas}\n').encode())
    print('Training file written!')
    data = df[df['test'] == True]
    with gzip.open(snakemake.output['test'], 'w') as f:
        f.write(b'FoldID\tEventID seq\tBound\t\n')
        for i, seq, clas in zip(data.index, data['seq'], data['class']):
            f.write((f'A\tseq_{i:07d}_peak\t{seq[:101]}\t{clas}\n').encode())
    print('Test file written!')

## preprocessing iDeepS
elif snakemake.params[0] == 'ideeps':
    print('Writing files for iDeepS...')
    data = df[df['test'] == False]
    with gzip.open(snakemake.output['train'], 'w') as csvfile:
        for descriptor, seq, clas in zip(data['descriptor'], data['seq'], data['class']):
            csvfile.write((descriptor + f'; class:{clas}\n').encode())
            csvfile.write((seq[:101] + '\n').encode())
    print('Training file written!')
    data = df[df['test'] == True]
    with gzip.open(snakemake.output['test'], 'w') as csvfile:
        for descriptor, seq, clas in zip(data['descriptor'], data['seq'], data['class']):
            csvfile.write((descriptor + f'; class:{clas}\n').encode())
            csvfile.write((seq[:101] + '\n').encode())
    print('Test file written!')

## preprocessing GraphProt2
elif snakemake.params[0] == 'graphprot2':
    print('Writing files for GraphProt2...')
    # needs two files: one with positives and one with negatives
    # no sequences with only N are allowed -> already removed
    # no duplicated chromosome names are allowed -> use index as identifier
    data = df[df['test'] == False]
    # write positive file
    with open(snakemake.output['positive'], 'w', newline='') as f:
        positives = data[data['class'] == 1]
        for descriptor, seq in zip(positives.index, positives['seq']):
            f.write(f'>{descriptor:05d}\n')
            f.write(seq + '\n')
    print('Positive File written!')

    # write negative file
    with open(snakemake.output['negative'], 'w', newline='') as f:
        negatives = data[data['class'] == 0]
        for descriptor, seq in zip(negatives.index, negatives['seq']):
            f.write('>' + str(descriptor) + '\n')
            f.write(seq + '\n')
    print('Negative file written!')

    # write test file
    data = df[df['test'] == True]
    with open(snakemake.output['test'], 'w', newline='') as f:
        for descriptor, seq in zip(data.index, data['seq']):
            f.write('>' + str(descriptor) + '\n')
            f.write(seq + '\n')

    print('Test file written!')
