import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import accuracy_score
from sklearn.metrics import average_precision_score

methods = snakemake.params['methods']

print(f'Creating report for {methods}...\n')

# load data
df = pd.read_csv(snakemake.input[0]).set_index(['fold', 'Sample ID'])
folds = snakemake.config['folds']

print('Data loaded')

### Calc auROC and accuracy ###

auROCs = dict([(method,[]) for method in methods])
accs = dict([(method,[]) for method in methods])
aps = dict([(method,[]) for method in methods])
for method in methods:
    for i in range(folds):

        # select data
        y_true = df.loc[i, 'class']
        predictions = df.loc[i, method]

        # convert predictions to labels 0 if < 0.5 else 1
        y_pred = predictions.mask(predictions < 0.5, 0)
        y_pred = y_pred.mask(y_pred != 0, 1)
        y_pred = y_pred.astype('int')

        acc_score = accuracy_score(y_true, y_pred)

        # calc ROC curve & AUC
        fpr, tpr, _ = roc_curve(y_true, predictions)
        roc_auc = auc(fpr, tpr)

        # calc average_precision
        average_precision = average_precision_score(y_true, predictions)

        accs[method].append(acc_score)
        auROCs[method].append(roc_auc)
        aps[method].append(average_precision)


### Create graphs ###
scale=snakemake.params['scale']
fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(scale*3,scale*1), sharey=True)

# Boxplot for ROC auc values over folds
ax1.boxplot(auROCs.values(),
            vert=True,
            patch_artist=True,
            labels=auROCs.keys())
ax1.set_title(f'auROC')

# Boxplot for accuracy values over folds
ax2.boxplot(accs.values(),
            vert=True,
            patch_artist=True,
            labels=accs.keys())
ax2.set_title(f'Accuracy')

# Boxplot for AP values over folds
ax3.boxplot(aps.values(),
            vert=True,
            patch_artist=True,
            labels=aps.keys())
ax3.set_title(f'Average Precision')

# define y axis
start = 0.4
ax1.set_ylim([start,1.])
ax1.set_yticks(np.arange(start,1.1, step=0.1))


fig.suptitle(f'{snakemake.wildcards.dataset} with {folds}-fold cross-validation')
fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=300)
