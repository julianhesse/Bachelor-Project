import pandas as pd
import numpy as np

print('Start preprocessing ...')
method = snakemake.params["method"]
prediction_range = snakemake.config['prediction_range'][method]
print(method + ':', prediction_range)

# load prediction values
if method == 'graphprot2':
    predictions = pd.read_csv(snakemake.input[0], header=None, sep='\t', index_col=0,
        squeeze=True)
    predictions.sort_index(inplace=True)
else:
    predictions = pd.read_csv(snakemake.input[0], header=None, sep='\t', squeeze=True)

print('Predictions loaded:')
print(predictions)

# transform to range (-1, 1)
if prediction_range != [-1, 1]:
    arr = predictions
    arr = arr - prediction_range[0]
    scale = prediction_range[1] - prediction_range[0]
    arr = (arr / scale) * 2 - 1
    
    predictions = arr

print('\nPredictions transformed:')
print(predictions)

# save transformed predictions
predictions.to_csv(snakemake.output[0])

print('\nSaved to:', snakemake.output[0])
print('\nDone!\n')
