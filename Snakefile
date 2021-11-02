import os

print(
"""
  ______                 _                          _    
  | ___ \               | |                        | |   
  | |_/ / ___ _ __   ___| |__  _ __ ___   __ _ _ __| | __
  | ___ \/ _ \ '_ \ / __| '_ \| '_ ` _ \ / _` | '__| |/ /
  | |_/ /  __/ | | | (__| | | | | | | | | (_| | |  |   < 
  \____/ \___|_| |_|\___|_| |_|_| |_| |_|\__,_|_|  |_|\_\\

"""
)

configfile: "config.yaml"

DATA_FILE="out/{dataset, [A-Za-z0-9_]+}/db.csv" # Database for data about the samples
DATA_FILE_="out/{dataset}/db.csv" # Database for data about the samples
DATASETS=os.listdir('datasets')
OUT_FOLDER=glob_wildcards('out/{dataset, [A-Za-z0-9_]+}').dataset
OUT_FOLDER_=glob_wildcards('out/{dataset}').dataset


#### preprocessing #####

rule preprocess_all:
    input:
        expand("out/{dataset}/data/{fold}_train_deepbind.seq.gz", dataset=["RBFOX2_HepG2_iDeepS"], fold=[i for i in range(5)]),
        expand("out/{dataset}/data/{fold}_test_deepbind.seq.gz", dataset=["RBFOX2_HepG2_iDeepS"], fold=[i for i in range(5)]),
        expand("out/{dataset}/data/{fold}_train_ideeps.fa.gz", dataset=["RBFOX2_HepG2_iDeepS"], fold=[i for i in range(5)]),
        expand("out/{dataset}/data/{fold}_test_ideeps.fa.gz", dataset=["RBFOX2_HepG2_iDeepS"], fold=[i for i in range(5)]),
        expand("out/{dataset}/data/{fold}_train_graphprot_positive.fasta", dataset=["RBFOX2_HepG2_iDeepS"], fold=[i for i in range(5)]),
        expand("out/{dataset}/data/{fold}_train_graphprot_negative.fasta", dataset=["RBFOX2_HepG2_iDeepS"], fold=[i for i in range(5)]),
        expand("out/{dataset}/data/{fold}_test_graphprot.fasta", dataset=["RBFOX2_HepG2_iDeepS"], fold=[i for i in range(5)])


rule preprocess_dataset_folds:
    input:
        positive=expand("datasets/{dataset}/{format}/fold-{fold}/positive.fold-{fold}.fasta",
                dataset="{dataset}", format="{format}",
                fold=[i for i in range(config['folds'])]),
        negative=expand("datasets/{dataset}/{format}/fold-{fold}/negative-1-2.fold-{fold}.fasta",
                dataset="{dataset}", format="{format}",
                fold=[i for i in range(config['folds'])])
    output:
        "out/{dataset, [A-Za-z0-9_]+}_{format, [A-Za-z]+}/db.csv"
    params:
        folds='manual'
    log:
        "logs/preprocess/{dataset}_{format}_to_csv.log"
    script:
        "scripts/preprocess_dataset.py"

rule preprocess_dataset:
    input:
        positive="datasets/{dataset}/positive.fasta",
        negative="datasets/{dataset}/negative-1-2.fasta"
    output:
        DATA_FILE
    params:
        folds=config['folds'],
        seed=462
    log:
        "logs/preprocess/{dataset}_to_csv.log"
    script:
        "scripts/preprocess_dataset.py"

rule preprocess_deepbind:
    input:
        data=DATA_FILE_,
    params:
        method="deepbind",
        fold="{wildcards.fold}"
    output:
        train=temp("out/{dataset}/data/{fold, [0-9]+}_train_deepbind.seq.gz"),
        test=temp("out/{dataset}/data/{fold, [0-9]+}_test_deepbind.seq.gz"),
    log:
        "logs/preprocess/{dataset}_fold-{fold}_preprocess_deepbind.log"
    script:
        "scripts/preprocess_for_methods.py"

rule preprocess_ideeps:
    input:
        data=DATA_FILE_,
    params:
        method="ideeps",
        fold="{wildcards.fold}"
    output:
        train=temp("out/{dataset}/data/{fold, [0-9]+}_train_ideeps.fa.gz"),
        test=temp("out/{dataset}/data/{fold, [0-9]+}_test_ideeps.fa.gz"),
    log:
        "logs/preprocess/{dataset}_fold-{fold}_preprocess_deepbind.log"
    script:
        "scripts/preprocess_for_methods.py"

rule preprocess_graphprot:
    # GraphProt and GraphProt2 can use the same format!
    # very nice :)
    input:
        data=DATA_FILE_,
    params:
        method="graphprot",
        fold="{wildcards.fold}"
    output:
        positive=temp("out/{dataset}/data/{fold, [0-9]+}_train_graphprot_positive.fasta"),
        negative=temp("out/{dataset}/data/{fold, [0-9]+}_train_graphprot_negative.fasta"),
        test=temp("out/{dataset}/data/{fold, [0-9]+}_test_graphprot.fasta")
    log:
        "logs/preprocess/{dataset}_fold-{fold}_preprocess_deepbind.log"
    script:
        "scripts/preprocess_for_methods.py"



#### train & predict #####


rule run_deepbind:
    input:
        train="out/{dataset}/data/{fold}_train_deepbind.seq.gz",
        test="out/{dataset}/data/{fold}_test_deepbind.seq.gz",
    params:
        "out/{dataset}/fold-{fold}/deepbind/"
    output:
        prediction="out/{dataset}/fold-{fold, [0-9]+}/deepbind/prediction.out"
    benchmark:
        "out/{dataset}/fold-{fold}/deepbind/benchmark.txt"
    threads: 3
    resources: gpu=1, mem_mb=4500, time_min=45, partition="gpu_p", qos="gpu"
    conda:
        "envs/deepbind.yaml"
    log:
        "logs/out/deepbind/{dataset}_fold-{fold, [0-9]+}_deepbind_run.log"
    script:
        "methods/DeepBind_with_Tensorflow/deepbind.py"

rule run_deepbind_pytorch:
    input:
        train="out/{dataset}/data/{fold}_train_deepbind.seq.gz",
        test="out/{dataset}/data/{fold}_test_deepbind.seq.gz"
    params:
        out="out/{dataset}/fold-{fold}/deepbind_pytorch/"
    output:
        prediction="out/{dataset}/fold-{fold, [0-9]+}/deepbind_pytorch/prediction.out"
    benchmark:
        "out/{dataset}/fold-{fold}/deepbind_pytorch/benchmark.txt"
    threads: 3
    resources: gpu=1, mem_mb=4500, time_min=45, partition="gpu_p", qos="gpu"
    conda:
        "envs/deepbind_pytorch.yaml"
    log:
        "logs/out/deepbind/{dataset}_fold-{fold, [0-9]+}_deepbind_run.log"
    shell:
        "python methods/deepbind_pytorch/deepbind.py {input.train} {input.test} {params.out} &> {log}"

rule run_ideeps:
    input:
        train="out/{dataset}/data/{fold}_train_ideeps.fa.gz",
        test="out/{dataset}/data/{fold}_test_ideeps.fa.gz",
    params:
        "out/{dataset}/fold-{fold}/ideeps/"
    output:
        model="out/{dataset}/fold-{fold, [0-9]+}/ideeps/model.pkl",
        prediction="out/{dataset}/fold-{fold, [0-9]+}/ideeps/prediction.out"
    benchmark:
        "out/{dataset}/fold-{fold}/ideeps/benchmark.txt"
    threads: 2
    resources: cpus=2, mem_mb=5000, time_min=120, gpu=0, partition="gpu_p", qos="gpu"
    conda:
        "envs/ideeps.yaml"
    log:
        "logs/out/ideeps/{dataset}_fold-{fold}_iDeepS_run.log"
    shell:
        "python scripts/run_ideeps.py {input.train} {input.test} {output.prediction} {params[0]} {wildcards.fold} &> {log}"

rule run_graphprot2:
    input:
        positive="out/{dataset}/data/{fold}_train_graphprot_positive.fasta",
        negative="out/{dataset}/data/{fold}_train_graphprot_negative.fasta",
        test="out/{dataset}/data/{fold}_test_graphprot.fasta"
    params:
        out="out/{dataset}/fold-{fold}/graphprot2",
        conainer="charliecloud"
        #conainer="singularity"
    output:
        model="out/{dataset}/fold-{fold, [0-9]+}/graphprot2/trained_model/final.model",
        prediction="out/{dataset}/fold-{fold, [0-9]+}/graphprot2/prediction/whole_site_scores.out"
    benchmark:
        "out/{dataset}/fold-{fold}/graphprot2/benchmark.txt"
    threads: 3
    resources: gpu=1, mem_mb=4000, time_min=60, partition="gpu_p", qos="gpu"
    conda:
        "envs/graphprot2.yaml"
    log:
        "logs/out/graphprot2/{dataset}_fold-{fold}_graphprot2_run.log"
    shell:
        """
        graphprot2 gt --in {input.positive} --neg-in {input.negative} --out {params.out}/train_data &> {log}
        graphprot2 train --in {params.out}/train_data --out {params.out}/trained_model &>> {log}

        graphprot2 gp --in {input.test} --train-in {params.out}/trained_model --out {params.out}/prediction_data &>> {log}
        graphprot2 predict --in {params.out}/prediction_data --train-in {params.out}/trained_model --out {params.out}/prediction --mode 1 &>> {log}
        """

rule run_graphprot:
    input:
        positive="out/{dataset}/data/{fold}_train_graphprot_positive.fasta",
        negative="out/{dataset}/data/{fold}_train_graphprot_negative.fasta",
        test="out/{dataset}/data/{fold}_test_graphprot.fasta"
    output:
        model="out/{dataset}/fold-{fold, [0-9]+}/graphprot/GraphProt.model",
        prediction="out/{dataset}/fold-{fold, [0-9]+}/graphprot/prediction.out"
    benchmark:
        "out/{dataset}/fold-{fold}/graphprot/benchmark.txt"
    threads: 3
    resources: mem_mb=3000, time_min=60, partition="interactive_gpu_p", qos="interactive_gpu"
    conda:
        "envs/graphprot.yaml"
    log:
        "logs/out/graphprot/{dataset}_fold-{fold}_graphprot_run.log"
    shell:
        """
        cd out/{wildcards.dataset}/fold-{wildcards.fold}/graphprot &> {log}
        GraphProt.pl --action train -fasta ../../../../{input.positive} \
            -negfasta ../../../../{input.negative} &> ../../../../{log}
        GraphProt.pl --action predict -model GraphProt.model -fasta ../../../../{input.test} &> ../../../../{log}
        mv GraphProt.predictions prediction.out &> ../../../../{log}
        """

rule output_graphtprot2:
    input:
        "out/{dataset}/fold-{fold}/graphprot2/prediction/whole_site_scores.out"
    output:
        temp("out/{dataset}/fold-{fold, [0-9]+}/graphprot2/prediction.out")
    shell:
        "cp {input} {output}"

#### preprocess, aggreagate predictions ####

rule preprocess_prediction:
    input:
        "out/{dataset}/fold-{fold}/{method}/prediction.out"
    params:
        method="{method}"
    output:
        temp("out/{dataset, [A-Za-z0-9_]+}/fold-{fold, [0-9]+}/results/{method}_prediction.out")
    script:
        "scripts/preprocess_predictions.py"

def prediction_files(wildcards):
    predictions = {}
    dataset = wildcards.dataset
    fold = wildcards.fold
    for method in config['methods']:
        predictions[method] = f"out/{dataset}/fold-{fold}/results/{method}_prediction.out"
    return predictions

rule aggregate_fold_predictions:
    input:
        unpack(prediction_files),
        data=DATA_FILE_,
        # deepbind="out/{dataset}/fold-{fold}/results/deepbind_prediction.out",
        # ideeps="out/{dataset}/fold-{fold}/results/ideeps_prediction.out",
        # graphprot="out/{dataset}/fold-{fold}/results/graphprot_prediction.out",
        # graphprot2="out/{dataset}/fold-{fold}/results/graphprot2_prediction.out",
    output:
        temp("out/{dataset, [A-Za-z0-9_]+}/fold-{fold, [0-9]+}/results/predictions.csv")
    script:
        "scripts/aggregate_predictions.py"

def all_fold_results(wildcards):
    """
    For every fold the prediction.csv file are saved in a dict according to the fold-number.
    """
    files = {}
    for i in range(config['folds']):
        prediction_csv = "out/{wildcards.dataset}/fold-{fold}/results/predictions.csv".format(wildcards=wildcards, fold=i)
        files[str(i)] = prediction_csv
    return files

rule aggregate_dataset_predictions:
    input:
        unpack(all_fold_results)
    output:
        "out/{dataset, [A-Za-z0-9_]+}/results/predictions.csv"
    script:
        "scripts/aggregate_all_predictions.py"


#### evaluate predictions ####

def all_cell_line_results(wildcards):
    """
    Collect all the `predictions.csv` files for one cell line.
    """
    files = {}
    rbps = glob_wildcards("out/{{rbp}}_{cell_line}_iDeepS/results/predictions.csv".format(
           cell_line=wildcards.cell_line
         )).rbp
    for rbp in rbps:
        prediction_csv = f"out/{rbp}_{wildcards.cell_line}_iDeepS/results/predictions.csv"
        files[rbp] = prediction_csv
    return files

rule postprocessing_roc_auc:
     input:
        # files=expand("out/{rbp}_{cell_line}_iDeepS/results/predictions.csv", rpb=RBPs, cell_line=wildcards.cell_line)
        unpack(all_cell_line_results)
     output:
        "out/results/{cell_line, [A-Za-z0-9]+}_{mode}.csv"
     script:
        "scripts/metric_dataframe.py"


rule classification_report:
    input:
        "out/{dataset}/results/predictions.csv"
    output:
        "out/{dataset, [A-Za-z0-9_]+}/reports/{method}_report.txt",
    params:
        method="{method}"
    script:
        "scripts/classification_report.py"

rule classification_graphs:
    input:
        "out/{dataset}/results/predictions.csv"
    output:
        "out/{dataset, [A-Za-z0-9_]+}/reports/{method}_roc_pr_curve.png"
    params:
        method="{method}",
        scale=5
    script:
        "scripts/classification_graphs.py"

rule performance_comp:
    input:
        "out/{dataset}/results/predictions.csv"
    output:
        "out/{dataset, [A-Za-z0-9_]+}/reports/performance_comp.png"
    params:
        methods=config['methods'],
        scale=4
    script:
        "scripts/performance_comp.py"



rule report_all:
    input:
        expand("out/{dataset}/reports/{method}_report.txt",
                dataset=config['datasets'],
                method=config['methods']),
        expand("out/{dataset}/reports/{method}_roc_pr_curve.png",
                dataset=config['datasets'],
                method=config['methods']),
        expand("out/{dataset}/reports/performance_comp.png",
                dataset=config['datasets'],
                method=config['methods'])

rule report_all_test:
    input:
        expand("out/Test/reports/{method}_report.txt",
                method=config['methods']),
        expand("out/Test/reports/{method}_roc_pr_curve.png",
                method=config['methods']),
        expand("out/Test/reports/performance_comp.png",
                method=config['methods'])

rule compress:
    input:
        "out/{dataset}/reports/performance_comp.png"
    output:
        "out/{dataset, [A-Za-z0-9_]+}.tar.gz"
    shell:
        """
        cd out
        tar -zcf {wildcards.dataset}.tar.gz {wildcards.dataset}
        """

rule compress_all:
    input:
        #expand("out/{dataset}.tar.gz", dataset=glob_wildcards(DATA_FILE).dataset)
        expand("out/{dataset}.tar.gz", dataset=config['datasets'])


#### Fancy graphs ####

rule scatter_plot:
    input:
        "out/results/{cell_line}_{mode}.csv"
    output:
        "out/plots/{cell_line, [A-Za-z]+[0-9]+}_{method_0}_{method_1}_scatterplot_{mode}.png"
    script:
        "scripts/scatter_plot.py"

rule box_plot:
     input:
        "out/results/{cell_line}_{mode}.csv"
     output:
        "out/plots/{cell_line}_boxplot_{mode}.png"
     script:
        "scripts/box_plot.py"

rule bar_plot:
     input:
        "out/results/{cell_line}_{mode}.csv"
     output:
        "out/plots/{cell_line}_barplot_{mode}.png"
     script:
        "scripts/bar_plot.py"

rule plots_roc_auc:
     input:
        expand("out/plots/{cell_line}_{method_0}_{method_1}_scatterplot_roc_auc.png",
            cell_line=config['cell_lines'],
            method_0=config['methods'],
            method_1=config['methods']
        ),
        expand("out/plots/{cell_line}_boxplot_roc_auc.png", cell_line=config['cell_lines'])

rule plots_ap:
     input:
        expand("out/plots/{cell_line}_{method_0}_{method_1}_scatterplot_ap.png",
            cell_line=config['cell_lines'],
            method_0=config['methods'],
            method_1=config['methods']
        ),
        expand("out/plots/{cell_line}_boxplot_ap.png", cell_line=config['cell_lines'])

def my_func(wildcards):
    print('Wildcards:')
    print(wildcards)
    for i in range(config['folds']):
        print("out/{wildcards.dataset}/fold-{fold, [0-9]+}/".format(wildcards=wildcards, fold=i))
    return {}

#import os

#DATASET_FOLDER="datasets/{dataset}"
#DATASETS=os.listdir('datasets')

#rule test:
    #input:
        #datasets=expand(DATASET_FOLDER, dataset=DATASETS)
        ##unpack(my_func),
        ##pred=expand('out/RBFOX2_HepG2_iDeepS/{method}/benchmark.txt', method=config['methods'])
    ##params:
        ##folds=glob_wildcards('out/{dataset}/fold-{fold, [0-9]+}/{method}/benchmark.txt').fold
    ##output:
        ##'out/{dataset}/done'
    #run:
        #print(input['datasets'])
        #print(list(input['datasets']))
        #print('Datasets:', DATASETS)
