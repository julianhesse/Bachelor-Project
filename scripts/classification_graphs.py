import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, roc_curve, auc
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay
from sklearn.metrics import precision_recall_curve, average_precision_score

method = snakemake.params['method']

print(f'Creating graphs for {method}...\n')

# load data
df = pd.read_csv(snakemake.input[0]).set_index(['fold', 'Sample ID'])
folds = snakemake.config['folds']

print('Data loaded')


### plot roc curve and precision recall curve ###

# create graph
# first row is for all graphs combined
# rows after that are for folds
scale=5
fig, axs = plt.subplots(folds + 1, 2, sharex='col', sharey='col', figsize=(scale*2, scale*(folds)))

# prepare plots
lw = 1.5

# roc plots
axs[0,0].plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')

axs[0,0].set_title(f'Receiver operating characteristic {method}')
axs[0,0].set_xlim([0.0, 1.0])
axs[0,0].set_ylim([0.0, 1.05])
axs[0,0].set_xlabel('False Positive Rate')
axs[0,0].set_ylabel('True Positive Rate')

# pr plots 
axs[0,1].set_title(f'2-class Precision-Recall curve {method}')
axs[0,1].set_xlim([-0.1, 1.1])
axs[0,1].set_ylim([0.0, 1.05])
axs[0,1].set_xlabel('Recall')
axs[0,1].set_ylabel('Precision')

for i in range(folds):
    
    # select data
    y_true = df.loc[i, 'class']
    predictions = df.loc[i, method]

    # calc ROC curve & AUC
    fpr, tpr, _ = roc_curve(y_true, predictions)
    roc_auc = auc(fpr, tpr)


    # calc PR curve & average_precision
    prec, recall, _ = precision_recall_curve(y_true, predictions)
    average_precision = average_precision_score(y_true, predictions)


    # plot ROC curve
    axs[0,0].plot(fpr, tpr, #color='darkorange',
             lw=lw, label='Fold %d' % i)
    axs[i+1,0].plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    axs[i+1,0].plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    axs[i+1,0].set_xlabel('False Positive Rate')
    axs[i+1,0].set_ylabel('True Positive Rate')
    axs[i+1,0].set_title(f'ROC Fold-{i}')
    axs[i+1,0].legend(loc="lower right")

    # plot PR curve
    axs[0,1].plot(recall, prec,
             lw=lw, label='Fold %d' % i)
    axs[i+1,1].plot(recall, prec,
             lw=lw, label='PR curve (AP = %0.2f)' % average_precision)
    axs[i+1,1].set_xlabel('Recall')
    axs[i+1,1].set_ylabel('Precision')
    axs[i+1,1].set_title(f'PR Fold-{i}')
    axs[i+1,1].legend(loc="lower left")


    #roc_display = RocCurveDisplay(fpr=fpr, tpr=tpr).plot(ax=ax3)
     
    #pr_display = PrecisionRecallDisplay(precision=prec, recall=recall).plot(ax=ax4)

# display legends
axs[0,0].legend(loc="lower right")
axs[0,1].legend(loc="lower left")

# fig.savefig(snakemake.output[0], dpi=300, bbox_inches='tight')
fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=300)

print(f'\nPlots saved to {snakemake.output[0]}')
print('\nDone!\n')
