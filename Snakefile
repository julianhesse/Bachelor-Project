configfile: "config.yaml"


#### preprocessing #####

rule preprocess_all:
    input:
        "preprocessed/RBFOX2_HepG2_iDeepS/train_deepbind.seq.gz",
        "preprocessed/RBFOX2_HepG2_iDeepS/test_deepbind.seq.gz",
        "preprocessed/RBFOX2_HepG2_iDeepS/train_ideeps.fa.gz",
        "preprocessed/RBFOX2_HepG2_iDeepS/test_ideeps.fa.gz",
        "preprocessed/RBFOX2_HepG2_iDeepS/train_graphprot2_positive.fasta",
        "preprocessed/RBFOX2_HepG2_iDeepS/train_graphprot2_negative.fasta",
        "preprocessed/RBFOX2_HepG2_iDeepS/test_graphprot2.fasta"

rule preprocess_dataset:
    input:
        positive="datasets/{dataset}/positive.fasta",
        negative="datasets/{dataset}/negative-1-2.fasta"
    output:
        "datasets/{dataset}/temp.csv"
    params:
        frac=0.25,
        seed=462
    log:
        "logs/preprocess/{dataset}_to_csv.log"
    script:
        "scripts/preprocess_dataset.py"

rule preprocess_deepbind:
    input:
        data="datasets/RBFOX2_HepG2_iDeepS/temp.csv",
    params:
        method="deepbind"
    output:
        train="preprocessed/RBFOX2_HepG2_iDeepS/train_deepbind.seq.gz",
        test="preprocessed/RBFOX2_HepG2_iDeepS/test_deepbind.seq.gz",
    script:
        "scripts/preprocess_for_methods.py"

rule preprocess_ideeps:
    input:
        data="datasets/RBFOX2_HepG2_iDeepS/temp.csv",
    params:
        method="ideeps"
    output:
        train="preprocessed/RBFOX2_HepG2_iDeepS/train_ideeps.fa.gz",
        test="preprocessed/RBFOX2_HepG2_iDeepS/test_ideeps.fa.gz",
    script:
        "scripts/preprocess_for_methods.py"

rule preprocess_graphprot2:
    input:
        data="datasets/RBFOX2_HepG2_iDeepS/temp.csv",
    params:
        method="graphprot2"
    output:
        positive="preprocessed/RBFOX2_HepG2_iDeepS/train_graphprot2_positive.fasta",
        negative="preprocessed/RBFOX2_HepG2_iDeepS/train_graphprot2_negative.fasta",
        test="preprocessed/RBFOX2_HepG2_iDeepS/test_graphprot2.fasta"
    script:
        "scripts/preprocess_for_methods.py"



#### train & predict #####

rule run_deepbind:
    input:
        train="preprocessed/RBFOX2_HepG2_iDeepS/train_deepbind.seq.gz",
        test="preprocessed/RBFOX2_HepG2_iDeepS/test_deepbind.seq.gz",
    params:
        "out/RBFOX2_HepG2_iDeepS/deepbind/"
    output:
        prediction="out/RBFOX2_HepG2_iDeepS/deepbind/prediction.out"
    conda:
        "envs/deepbind.yaml"
    log:
        "logs/out/deepbind/RBFOX2_HepG2_iDeepS_run.log"
    script:
        "methods/DeepBind_with_Tensorflow/deepbind.py"


rule run_ideeps:
    input:
        train="preprocessed/RBFOX2_HepG2_iDeepS/train_ideeps.fa.gz",
        test="preprocessed/RBFOX2_HepG2_iDeepS/test_ideeps.fa.gz",
    params:
        "out/RBFOX2_HepG2_iDeepS/ideeps/"
    output:
        model="out/RBFOX2_HepG2_iDeepS/ideeps/model.pkl",
        prediction="out/RBFOX2_HepG2_iDeepS/ideeps/prediction.out"
    conda:
        "envs/ideeps.yaml"
    log:
        "logs/out/ideeps/RBFOX2_HepG2_iDeepS_run.log"
    script:
        "scripts/run_ideeps.py"

rule run_graphprot2:
    input:
        positive="preprocessed/RBFOX2_HepG2_iDeepS/train_graphprot2_positive.fasta",
        negative="preprocessed/RBFOX2_HepG2_iDeepS/train_graphprot2_negative.fasta",
        test="preprocessed/RBFOX2_HepG2_iDeepS/test_graphprot2.fasta"
    params:
        "out/RBFOX2_HepG2_iDeepS/graphprot2"
    output:
        #model="out/RBFOX2_HepG2_iDeepS/graphprot2/trained_model/final.model",
        prediction="out/RBFOX2_HepG2_iDeepS/graphprot2/prediction/whole_site_scores.out"
    # conda:
    #     "envs/graphprot2.yaml"
    log:
        "logs/out/graphprot2/RBFOX2_HepG2_iDeepS_run.log"
    script:
        "scripts/run_graphprot2.py"

#### aggreagate predictions & evaluate ####

rule aggregate_predictions:
    input:
        dataset="datasets/{dataset}/temp.csv",
        deepbind="out/{dataset}/deepbind/prediction.out",
        ideeps="out/{dataset}/ideeps/prediction.out",
        graphprot2="out/{dataset}/graphprot2/prediction/whole_site_scores.out"
    output:
        "results/{dataset}/test_results.csv"
    script:
        "scripts/aggregate_predictions.py"
