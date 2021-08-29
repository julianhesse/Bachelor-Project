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
        "out/{dataset}/temp.csv"
    params:
        frac=0.25,
        seed=462
    log:
        "logs/preprocess/{dataset}_to_csv.log"
    script:
        "scripts/preprocess_dataset.py"

rule preprocess_deepbind:
    input:
        data="out/{dataset}/temp.csv",
    params:
        method="deepbind"
    output:
        train="out/{dataset}/data/train_deepbind.seq.gz",
        test="out/{dataset}/data/test_deepbind.seq.gz",
    script:
        "scripts/preprocess_for_methods.py"

rule preprocess_ideeps:
    input:
        data="out/{dataset}/temp.csv",
    params:
        method="ideeps"
    output:
        train="out/{dataset}/data/train_ideeps.fa.gz",
        test="out/{dataset}/data/test_ideeps.fa.gz",
    script:
        "scripts/preprocess_for_methods.py"

rule preprocess_graphprot:
    # GraphProt and GraphProt2 can use the same format!
    # very nice :)
    input:
        data="out/{dataset}/temp.csv",
    params:
        method="graphprot"
    output:
        positive="out/{dataset}/data/train_graphprot_positive.fasta",
        negative="out/{dataset}/data/train_graphprot_negative.fasta",
        test="out/{dataset}/data/test_graphprot.fasta"
    script:
        "scripts/preprocess_for_methods.py"



#### train & predict #####

rule run_deepbind:
    input:
        train="out/{dataset}/data/train_deepbind.seq.gz",
        test="out/{dataset}/data/test_deepbind.seq.gz",
    params:
        "out/{dataset}/deepbind/"
    output:
        prediction="out/{dataset}/deepbind/prediction.out"
    benchmark:
        "out/{dataset}/deepbind/benchmark.txt"
    conda:
        "envs/deepbind.yaml"
    log:
        "logs/out/deepbind/{dataset}_deepbind_run.log"
    script:
        "methods/DeepBind_with_Tensorflow/deepbind.py"


rule run_ideeps:
    input:
        train="out/{dataset}/data/train_ideeps.fa.gz",
        test="out/{dataset}/data/test_ideeps.fa.gz",
    params:
        "out/{dataset}/ideeps/"
    output:
        model="out/{dataset}/ideeps/model.pkl",
        prediction="out/{dataset}/ideeps/prediction.out"
    benchmark:
        "out/{dataset}/ideeps/benchmark.txt"
    conda:
        "envs/ideeps.yaml"
    log:
        "logs/out/ideeps/{dataset}_iDeepS_run.log"
    script:
        "scripts/run_ideeps.py"

rule run_graphprot2:
    input:
        positive="out/{dataset}/data/train_graphprot_positive.fasta",
        negative="out/{dataset}/data/train_graphprot_negative.fasta",
        test="out/{dataset}/data/test_graphprot.fasta"
    params:
        "out/{dataset}/graphprot2",
        conainer="charliecloud"
        #conainer="singularity"
    output:
        model="out/{dataset}/graphprot2/trained_model/final.model",
        prediction="out/{dataset}/graphprot2/prediction/whole_site_scores.out"
    benchmark:
        "out/{dataset}/graphprot2/benchmark.txt"
    log:
        "logs/out/graphprot2/{dataset}_graphprot2_run.log"
    script:
        "scripts/run_graphprot2.py"

rule run_graphprot:
    input:
        positive="out/{dataset}/data/train_graphprot_positive.fasta",
        negative="out/{dataset}/data/train_graphprot_negative.fasta",
        test="out/{dataset}/data/test_graphprot.fasta"
    output:
        model="out/{dataset}/graphprot/GraphProt.model",
        prediction="out/{dataset}/graphprot/prediction.out"
    benchmark:
        "out/{dataset}/graphprot/benchmark.txt"
    conda:
        "envs/graphprot.yaml"
    log:
        "logs/out/graphprot2/{dataset}_graphprot_run.log"
    shell:
        """
        cd out/{wildcards.dataset}/graphprot
        GraphProt.pl --action train -fasta ../../../{input.positive} \
            -negfasta ../../../{input.negative}
        GraphProt.pl --action predict -model GraphProt.model -fasta ../../../{input.test}
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
