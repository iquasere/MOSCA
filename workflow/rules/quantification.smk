rule quantification:
    input:
        expand("{output}/Preprocess/Trimmomatic/quality_trimmed_{name}{fr}.fq", output=OUTPUT, name=not_mp_exps["Name"],
            fr=(['_forward_paired', '_reverse_paired'] if EXPS["Files"].str.contains(',').tolist() else '')),
        expand("{output}/Assembly/{sample}/contigs.fasta", output=OUTPUT, sample=set(EXPS["Sample"])) if config['do_assembly'] else [],
        expand("{output}/Annotation/{sample}/fgs.ffn", output=OUTPUT, sample=set(EXPS["Sample"]))
    output:
        expand("{output}/Quantification/{name}.readcounts", output=OUTPUT, name=set(not_mp_exps['Name'])),
        expand("{output}/Quantification/{sample}_mg.readcounts", output=OUTPUT, sample=set(mg_exps['Sample'])) if len(mg_exps) > 0 else [],
        expand("{output}/Quantification/{sample}_mg_norm.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])) if (len(mg_exps) > 0 or not config['assembly']) else [],
        expand("{output}/Quantification/{sample}_mt.readcounts", output=OUTPUT, sample=set(mt_exps['Sample'])) if len(mt_exps) > 0 else [],
        expand("{output}/Quantification/{sample}_mt_norm.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])) if (len(mt_exps) > 0 or not config['assembly']) else []
    threads:
        config["threads"]
    params:
        output = OUTPUT,
        exps = f"{OUTPUT}/exps.tsv",
        did_assembly = config['do_assembly']
    conda:
        "../envs/quantification.yaml"
    script:
        "../scripts/quantification.py"
