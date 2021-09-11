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
    resources: gpu=1, mem_mb=4500, time_min=40
    conda:
        "envs/deepbind.yaml"
    log:
        "logs/out/deepbind/{dataset}_fold-{fold, [0-9]+}_deepbind_run.log"
    script:
        "methods/DeepBind_with_Tensorflow/deepbind.py"

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
    threads: 3
    resources: mem_mb=10000, time_min=60
    conda:
        "envs/ideeps.yaml"
    log:
        "logs/out/ideeps/{dataset}_fold-{fold}_iDeepS_run.log"
    script:
        "scripts/run_ideeps.py"

rule run_graphprot2:
    input:
        positive="out/{dataset}/data/{fold}_train_graphprot_positive.fasta",
        negative="out/{dataset}/data/{fold}_train_graphprot_negative.fasta",
        test="out/{dataset}/data/{fold}_test_graphprot.fasta"
    params:
        "out/{dataset}/fold-{fold}/graphprot2",
        conainer="charliecloud"
        #conainer="singularity"
    output:
        model="out/{dataset}/fold-{fold, [0-9]+}/graphprot2/trained_model/final.model",
        prediction="out/{dataset}/fold-{fold, [0-9]+}/graphprot2/prediction/whole_site_scores.out"
    benchmark:
        "out/{dataset}/fold-{fold}/graphprot2/benchmark.txt"
    threads: 3
    resources: mem_mb=3000, time_min=300
    log:
        "logs/out/graphprot2/{dataset}_fold-{fold}_graphprot2_run.log"
    script:
        "scripts/run_graphprot2.py"

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
    resources: mem_mb=8000, time_min=60
    conda:
        "envs/graphprot.yaml"
    log:
        "logs/out/graphprot/{dataset}_fold-{fold}_graphprot_run.log"
    shell:
        """
        cd out/{wildcards.dataset}/fold-{wildcards.fold}/graphprot
        GraphProt.pl --action train -fasta ../../../../{input.positive} \
            -negfasta ../../../../{input.negative}
        GraphProt.pl --action predict -model GraphProt.model -fasta ../../../../{input.test}
        mv GraphProt.predictions prediction.out
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
        temp("out/{dataset}/fold-{fold, [0-9]+}/results/{method}_prediction.out")
    script:
        "scripts/preprocess_predictions.py"


rule aggregate_fold_predictions:
    input:
        data=DATA_FILE_,
        deepbind="out/{dataset}/fold-{fold}/results/deepbind_prediction.out",
        ideeps="out/{dataset}/fold-{fold}/results/ideeps_prediction.out",
        graphprot="out/{dataset}/fold-{fold}/results/graphprot_prediction.out",
        graphprot2="out/{dataset}/fold-{fold}/results/graphprot2_prediction.out"
    output:
        temp("out/{dataset}/fold-{fold, [0-9]+}/results/predictions.csv")
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
        "out/{dataset}/results/predictions.csv"
    script:
        "scripts/aggregate_all_predictions.py"


#### evaluate predictions ####


rule classification_report:
    input:
        "out/{dataset}/results/predictions.csv"
    output:
        "out/{dataset}/reports/{method}_report.txt",
    params:
        method="{method}"
    script:
        "scripts/classification_report.py"

rule classification_graphs:
    input:
        "out/{dataset}/results/predictions.csv"
    output:
        "out/{dataset}/reports/{method}_roc_pr_curve.png"
    params:
        method="{method}",
        scale=5
    script:
        "scripts/classification_graphs.py"

rule performance_comp:
    input:
        "out/{dataset}/results/predictions.csv"
    output:
        "out/{dataset}/reports/performance_comp.png"
    params:
        methods=config['methods'],
        scale=4
    script:
        "scripts/performance_comp.py"


rule report_all:
    input:
        expand("out/RBFOX2_HepG2_iDeepS/reports/{method}_report.txt",
                method=config['methods']),
        expand("out/RBFOX2_HepG2_iDeepS/reports/{method}_roc_pr_curve.png",
                method=config['methods']),
        expand("out/RBFOX2_HepG2_iDeepS/reports/performance_comp.png",
                method=config['methods'])

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
