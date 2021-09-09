import pandas as pd
import gzip
import logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
        format='%(asctime)s %(name)s - %(levelname)s:%(message)s')

df = pd.read_csv(snakemake.input[0])

def write_deepbind(data, file):
    with gzip.open(file, 'w') as f:
        f.write(b'FoldID\tEventID seq\tBound\t\n')
        for i, seq, clas in zip(data.index, data['seq'], data['class']):
            f.write((f'A\tseq_{i:07d}_peak\t{seq[:101]}\t{clas}\n').encode())
    
def write_ideeps(data, file):
    with gzip.open(file, 'w') as csvfile:
        for descriptor, seq, clas in zip(data['descriptor'], data['seq'], data['class']):
            csvfile.write((descriptor + f'; class:{clas}\n').encode())
            csvfile.write((seq[:101] + '\n').encode())

def write_graphprot(data, file):
    with open(file, 'w', newline='') as f:
        for descriptor, seq in zip(data.index, data['seq']):
            f.write(f'>{descriptor:05d}\n')
            f.write(seq + '\n')

## preprocessing DeepBind
if snakemake.params['method'] == 'deepbind':
    logging.info('Writing files for DeepBind...')

    data = df[df['fold'] != int(snakemake.wildcards['fold'])]
    write_deepbind(data, snakemake.output['train'])
    logging.info('Training file written!')

    data = df[df['fold'] == int(snakemake.wildcards['fold'])]
    write_deepbind(data, snakemake.output['test'])
    logging.info('Test file written!')

## preprocessing iDeepS
elif snakemake.params['method'] == 'ideeps':
    logging.info('Writing files for iDeepS...')

    data = df[df['fold'] != int(snakemake.wildcards['fold'])]
    write_ideeps(data, snakemake.output['train'])
    logging.info('Training file written!')

    data = df[df['fold'] == int(snakemake.wildcards['fold'])]
    print(snakemake.wildcards['fold'], type(snakemake.wildcards['fold']))
    print(data)
    print(df)
    write_ideeps(data, snakemake.output['test'])
    logging.info('Test file written!')

## preprocessing GraphProt2
elif snakemake.params['method'] == 'graphprot':
    logging.info('Writing files for GraphProt and GraphProt2...')
    # needs two files: one with positives and one with negatives
    # no sequences with only N are allowed -> already removed
    # no duplicated chromosome names are allowed -> use index as identifier

    # write positive file
    data = df[df['fold'] != int(snakemake.wildcards['fold'])]
    data = data[data['class'] == 1]
    write_graphprot(data, snakemake.output['positive'])
    logging.info('Positive File written!')

    # write negative file
    data = df[df['fold'] != int(snakemake.wildcards['fold'])]
    data = data[data['class'] == 0]
    write_graphprot(data, snakemake.output['negative'])
    logging.info('Negative file written!')

    # write test file
    data = df[df['fold'] == int(snakemake.wildcards['fold'])]
    write_graphprot(data, snakemake.output['test'])
    with open(snakemake.output['test'], 'w', newline='') as f:
        for descriptor, seq in zip(data.index, data['seq']):
            f.write('>' + str(descriptor) + '\n')
            f.write(seq + '\n')

    logging.info('Test file written!')
