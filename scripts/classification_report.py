import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, roc_curve, auc
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay
from sklearn.metrics import precision_recall_curve, average_precision_score

method = snakemake.params['method']

print('Creating report for', method, '...\n')

# load data
df = pd.read_csv(snakemake.input[0])

print('Data loaded')

# create classification report for given method
y_true = df['class']
predictions = df[method]

## convert predictions to labels 0 if < 0.5 else 1
y_pred = predictions.mask(predictions < 0.5, 0)
y_pred = y_pred.mask(y_pred >= 0.5, 1)
y_pred = y_pred.astype('int')

# y_true = y_true.astype('float')

print(y_true)
print(y_pred)
print(predictions)

## create report with scit-learn
report = classification_report(y_true, y_pred, labels=[0,1])
print('Report created:')
print(report)

with open(snakemake.output[0], 'w') as f:
    f.write(report)

## plot roc curve and precision recall curve
# def true_false_positive(threshold_vector, y_test):
#     true_positive = np.equal(threshold_vector, 1) & np.equal(y_test, 1)
#     true_negative = np.equal(threshold_vector, 0) & np.equal(y_test, 0)
#     false_positive = np.equal(threshold_vector, 1) & np.equal(y_test, 0)
#     false_negative = np.equal(threshold_vector, 0) & np.equal(y_test, 1)
# 
#     tpr = true_positive.sum() / (true_positive.sum() + false_negative.sum())
#     fpr = false_positive.sum() / (false_positive.sum() + true_negative.sum())
# 
#     return tpr, fpr
# 
# def roc_from_scratch(probabilities, y_test, partitions=100):
#     roc = np.array([])
#     for i in range(partitions + 1):
#         
#         threshold_vector = np.greater_equal(probabilities, i / partitions).astype(int)
#         tpr, fpr = true_false_positive(threshold_vector, y_test)
#         roc = np.append(roc, [fpr, tpr])
#         
#     return roc.reshape(-1, 2)

# roc = roc_curve(y_true, y_pred)
# fpr = ROC[:,0]
# tpr = ROC[:,1]
fpr, tpr, _ = roc_curve(y_true, predictions)
roc_auc = auc(fpr, tpr)


prec, recall, _ = precision_recall_curve(y_true, predictions)
average_precision = average_precision_score(y_true, predictions)

fig, ax = plt.subplots(2, 2, figsize=(16, 16))

ax1, ax2, ax3, ax4 = ax.flat

# plt.figure()
lw = 1.5
ax1.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
ax1.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
ax1.set_xlim([0.0, 1.0])
ax1.set_ylim([0.0, 1.05])
ax1.set_xlabel('False Positive Rate')
ax1.set_ylabel('True Positive Rate')
ax1.set_title(f'Receiver operating characteristic {method}')
ax1.legend(loc="lower right")
#plt.show()

ax2.plot(recall, prec,
         lw=lw, label='PR curve (AP = %0.2f)' % average_precision)
ax2.set_xlim([-0.1, 1.1])
# ax2.set_ylim([0.0, 1.05])
ax2.set_xlabel('Recall')
ax2.set_ylabel('Precision')
ax2.set_title(f'2-class Precision-Recall curve {method}')
ax2.legend(loc="lower left")


roc_display = RocCurveDisplay(fpr=fpr, tpr=tpr).plot(ax=ax3)
 
pr_display = PrecisionRecallDisplay(precision=prec, recall=recall).plot(ax=ax4)

fig.savefig(snakemake.output[1], dpi=300)


print('\nPlots created')
print('\nDone!\n')
