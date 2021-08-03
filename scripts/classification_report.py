import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, roc_curve, auc

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

## create report with scit-learn
report = classification_report(y_true, y_pred, labels=[0,1])
print('Report created:')
print(report)

with open(snakemake.output[0], 'w') as f:
    f.write(report)

## plot roc curve
fpr, tpr, _ = roc_curve(y_true, y_pred)
roc_auc = auc(fpr, tpr)

plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(f'Receiver operating characteristic {method}')
plt.legend(loc="lower right")
#plt.show()
plt.savefig(snakemake.output[1])

print('\nPlots created')
print('\nDone!\n')
