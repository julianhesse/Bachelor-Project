import pandas as pd
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import average_precision_score

rbp_files = dict(snakemake.input)

folds = snakemake.config['folds']

# mode can be either 'roc_auc' or 'ap'
mode = snakemake.wildcards['mode']

methods = snakemake.config['methods']

## calc auc of roc-curve ##
def calc_roc_auc(y_true, predictions):
    fpr, tpr, _ = roc_curve(y_true, predictions)
    roc_auc = auc(fpr, tpr)

    return roc_auc

## calc average precision ##
def calc_ap(y_true, predictions):
    average_precision = average_precision_score(y_true, predictions)

    return average_precision



## iterate over all files ##
file_values = {}

for rbp, f_name in rbp_files.items():
    df = pd.read_csv(f_name).set_index(['fold', 'Sample ID'])

    values = dict([(method, []) for method in methods])

    for method in methods:
        for i in range(folds):
            # select data
            y_true = df.loc[i, 'class']
            predictions = df.loc[i, method]

            # calc values
            if mode == 'roc_auc':
                values[method].append(calc_roc_auc(y_true, predictions))
            elif mode == 'ap':
                values[method].append(calc_ap(y_true, predictions))
            else:
                raise Exception("No mode specified!!")

    df_rbp = pd.DataFrame(values)

    print(df_rbp)

    file_values[rbp] = df_rbp


## create multiindex dataframe ##
# The wanted dataframe:
# RBP Fold deepbind ideeps graphprot
# rbp    0    value  value     value
#        1    value  value     value
metric_df = pd.concat(file_values, axis=0, names=['RBP', 'fold'])

## add mean and std ##
metric_df = pd.concat(
    [metric_df, metric_df.groupby(level=0).agg(['mean', 'std']).stack(1)]
)


print(metric_df)
metric_df.to_csv(snakemake.output[0])
