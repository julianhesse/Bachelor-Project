configfile: "config.yaml"


#### preprocessing #####

DB="out/{dataset}/db.csv" # Database for data about the samples

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
        DB
    params:
        folds=config['folds'],
        seed=462
    log:
        "logs/preprocess/{dataset}_to_csv.log"
    script:
        "scripts/preprocess_dataset.py"

rule preprocess_deepbind:
    input:
        data=DB,
    params:
        method="deepbind",
        fold="{wildcards.fold}"
    output:
        train=temp("out/{dataset}/data/{fold}_train_deepbind.seq.gz"),
        test=temp("out/{dataset}/data/{fold}_test_deepbind.seq.gz"),
    log:
        "logs/preprocess/{dataset}_fold-{fold}_preprocess_deepbind.log"
    script:
        "scripts/preprocess_for_methods.py"

rule preprocess_ideeps:
    input:
        data=DB,
    params:
        method="ideeps",
        fold="{wildcards.fold}"
    output:
        train=temp("out/{dataset}/data/{fold}_train_ideeps.fa.gz"),
        test=temp("out/{dataset}/data/{fold}_test_ideeps.fa.gz"),
    log:
        "logs/preprocess/{dataset}_fold-{fold}_preprocess_deepbind.log"
    script:
        "scripts/preprocess_for_methods.py"

rule preprocess_graphprot:
    # GraphProt and GraphProt2 can use the same format!
    # very nice :)
    input:
        data=DB,
    params:
        method="graphprot",
        fold="{wildcards.fold}"
    output:
        positive=temp("out/{dataset}/data/{fold}_train_graphprot_positive.fasta"),
        negative=temp("out/{dataset}/data/{fold}_train_graphprot_negative.fasta"),
        test=temp("out/{dataset}/data/{fold}_test_graphprot.fasta")
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
        prediction="out/{dataset}/fold-{fold}/deepbind/prediction.out"
    benchmark:
        "out/{dataset}/fold-{fold}/deepbind/benchmark.txt"
    conda:
        "envs/deepbind.yaml"
    log:
        "logs/out/deepbind/{dataset}_fold-{fold}_deepbind_run.log"
    script:
        "methods/DeepBind_with_Tensorflow/deepbind.py"


rule run_ideeps:
    input:
        train="out/{dataset}/data/{fold}_train_ideeps.fa.gz",
        test="out/{dataset}/data/{fold}_test_ideeps.fa.gz",
    params:
        "out/{dataset}/fold-{fold}/ideeps/"
    output:
        model="out/{dataset}/fold-{fold}/ideeps/model.pkl",
        prediction="out/{dataset}/fold-{fold}/ideeps/prediction.out"
    benchmark:
        "out/{dataset}/fold-{fold}/ideeps/benchmark.txt"
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
        model="out/{dataset}/fold-{fold}/graphprot2/trained_model/final.model",
        prediction="out/{dataset}/fold-{fold}/graphprot2/prediction/whole_site_scores.out"
    benchmark:
        "out/{dataset}/fold-{fold}/graphprot2/benchmark.txt"
    log:
        "logs/out/graphprot2/{dataset}_fold-{fold}_graphprot2_run.log"
    script:
        "scripts/run_graphprot2.py"

rule run_graphprot:
    input:
        positive="out/{dataset}/data/train_graphprot_positive.fasta",
        negative="out/{dataset}/data/train_graphprot_negative.fasta",
        test="out/{dataset}/data/test_graphprot.fasta"
    output:
        model="out/{dataset}/fold-{fold}/graphprot/GraphProt.model",
        prediction="out/{dataset}/fold-{fold}/graphprot/prediction.out"
    benchmark:
        "out/{dataset}/fold-{fold}/graphprot/benchmark.txt"
    conda:
        "envs/graphprot.yaml"
    log:
        "logs/out/graphprot2/{dataset}_fold-{fold}_graphprot_run.log"
    shell:
        """
        cd out/{wildcards.dataset}/graphprot/fold-{wildcards.fold}
        GraphProt.pl --action train -fasta ../../../../{input.positive} \
            -negfasta ../../../../{input.negative}
        GraphProt.pl --action predict -model GraphProt.model -fasta ../../../../{input.test}
        mv GraphProt.predictions prediction.out
        """

rule output_graphtprot2:
    input:
        "out/{dataset}/graphprot2/prediction/whole_site_scores.out"
    output:
        temp("out/{dataset}/graphprot2/prediction.out")
    shell:
        "cp {input} {output}"

#### preprocess, aggreagate predictions ####

rule preprocess_prediction:
    input:
        "out/{dataset}/{method}/prediction.out"
    params:
        method="{method}"
    output:
        "out/{dataset}/results/{method}_prediction.out"
    script:
        "scripts/preprocess_predictions.py"

rule aggregate_predictions:
    input:
        dataset="out/{dataset}/temp.csv",
        deepbind="out/{dataset}/results/deepbind_prediction.out",
        ideeps="out/{dataset}/results/ideeps_prediction.out",
        graphprot="out/{dataset}/results/graphprot_prediction.out",
        graphprot2="out/{dataset}/results/graphprot2_prediction.out"
    output:
        "out/{dataset}/results/results.csv"
    script:
        "scripts/aggregate_predictions.py"

#### evaluate predictions ####
rule classification_report:
    input:
        "out/{dataset}/results/results.csv"
    output:
        "out/{dataset}/reports/{method}_report.txt",
        "out/{dataset}/reports/{method}_roc_pr_curve.png"
    params:
        method="{method}"
    script:
        "scripts/classification_report.py"

rule report_all:
    input:
        expand("out/RBFOX2_HepG2_iDeepS/reports/{method}_report.txt", method=config['methods'])
