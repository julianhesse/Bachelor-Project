import pandas as pd

## get the data from dataset
# but only the test set
df = pd.read_csv(snakemake.input['dataset'], index_col=0)

df = df[df['test'] == True]

## get the prediction values
methods = snakemake.config['methods']

print('Aggregation of predictions ...')
for method in methods:
    print('Loading predictions of', method, ':')
    print(snakemake.input[method], '\n')
    prediction = pd.read_csv(snakemake.input[method], squeeze=True)

    print(prediction)
    prediction.index = df.index
    print(prediction)
    print(df)
    df[method] = prediction

print(df)
print('\nAggreagation finished!\n')


## save results to csv for further analysis
df.to_csv(snakemake.output[0])
print('Written data to csv at:', snakemake.output[0], '\n')
