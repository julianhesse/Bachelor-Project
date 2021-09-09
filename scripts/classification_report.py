import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, roc_curve, auc
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay
from sklearn.metrics import precision_recall_curve, average_precision_score

method = snakemake.params['method']

print(f'Creating report for {method}...\n')

# load data
df = pd.read_csv(snakemake.input[0]).set_index(['fold', 'Sample ID'])
folds = snakemake.config['folds']

print('Data loaded')

### classification reports ###

reports = []

for i in range(folds):

    # select data
    y_true = df.loc[i, 'class']
    predictions = df.loc[i, method]

    # convert predictions to labels 0 if < 0.5 else 1
    y_pred = predictions.mask(predictions < 0.5, 0)
    y_pred = y_pred.mask(y_pred != 0, 1)
    y_pred = y_pred.astype('int')

    # y_true = y_true.astype('float')

    print(y_true)
    print(y_pred)
    print(predictions)

    # create report with scit-learn
    report = classification_report(y_true, y_pred, labels=[0,1])
    reports.append(report)
    print('Report created:')
    print(report)
    print()

with open(snakemake.output[0], 'w') as f:
    first = True
    for i, r in enumerate(reports):
        if first:
            first = False
        else:
            f.write("\n")
        f.write(f'Fold {i}:\n')
        f.write(r)

print(f'Reports written to: {snakemake.output[0]}')
print('\nDone!\n')
