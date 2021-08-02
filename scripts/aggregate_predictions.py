import pandas as pd

## get the data from dataset
# but only the test set
df = pd.read_csv(snakemake.input['dataset'], index_col=0)

df = df[df['test'] == True]

## get the prediction values
models = snakemake.config['models']

print('Aggregation of predictions ...')
for model in snakemake.config['models']:
    print('Loading predictions of', model, ':')
    print(snakemake.input[model], '\n')
    if model == 'graphprot2':
        prediction = pd.read_csv(snakemake.input[model], header=None, sep='\t', index_col=0,
        squeeze=True)
    else:
        prediction = pd.read_csv(snakemake.input[model], header=None, sep='\t', squeeze=True)

    prediction.index = df.index
    df[model] = prediction

print(df)
print('\nAggreagation finished!\n')


## save results to csv for further analysis
df.to_csv(snakemake.output[0])
print('Written data to csv at:', snakemake.output[0], '\n')
