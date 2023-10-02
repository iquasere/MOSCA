rule make_dea_input:
    input:
        expand("{output}/Quantification/mt.readcounts", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/mp_normalized.tsv", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        f"{OUTPUT}/DE_analysis/dea_input.tsv"
    threads:
        1
    run:
        result = pd.DataFrame(columns=['sseqid'])
        for file in input:
            df = pd.read_csv(file, sep='\t')
            sample = os.path.basename(os.path.dirname(file))
            # parse only first two columns of blast
            blast = pd.read_csv(
                f"{OUTPUT}/Annotation/{sample}/aligned.blast", sep='\t', usecols=[0,1], names=['qseqid', 'sseqid'])
            df = pd.merge(df, blast, how="left", on="qseqid")
            df = df.groupby('sseqid')[df.columns.tolist()[1:-1]].sum()
            # remove rows of df with less than 2 non-missing values
            df = df[(df > 0.0).sum(axis=1) > 1]
            if len(result) > 0:
                # inner join to keep all samples from having more than 2 non-missing values
                result = pd.merge(result, df, how="inner", on="sseqid")
            else:
                result = pd.merge(result, df, how="outer", on="sseqid")
        result.replace(0.0, np.nan).to_csv(output[0], sep='\t', index=False)

rule de_analysis:
    input:
        f"{OUTPUT}/DE_analysis/dea_input.tsv"
    output:
        f"{OUTPUT}/DE_analysis/condition_treated_results.tsv",
    threads:
        1
    params:
        conditions = (lambda wildcards, input: ",".join(not_mg_exps['Condition'].tolist())),
        foldchange = config["minimum_differential_expression"],
        fdr = config["false_discovery_rate"],
        output = (lambda wildcards, input: os.path.dirname(input[0])),
        datatype = (lambda wildcards, input: "rna_seq" if len(set(mt_exps['Sample'])) > 0 else "proteomics")
    conda:
        "../envs/de_analysis.yaml"
    script:
        "../scripts/de_analysis.R"
