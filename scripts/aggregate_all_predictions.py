import pandas as pd

print(snakemake.input)
fold_prediction_csvs = dict(snakemake.input)
print(fold_prediction_csvs)
folds = int(snakemake.config['folds'])

frames = []
for i in range(folds):
    df = pd.read_csv(fold_prediction_csvs[str(i)], index_col=0)
    frames.append(df)
    print(fold_prediction_csvs[str(i)])
    print(df)


all_predictions = pd.concat(frames, axis=0,
        keys=[i for i in range(folds)], names=['fold', 'Sample ID'])

print(all_predictions)
all_predictions.to_csv(snakemake.output[0])
