rule quantification:
    input:
        expand("{output}/Preprocess/Trimmomatic/quality_trimmed_{name}{fr}.fq", output=OUTPUT, name=not_mp_exps["Name"],
            fr=(['_forward_paired', '_reverse_paired'] if EXPS["Files"].str.contains(',').tolist() else '')),
        expand("{output}/Assembly/{sample}/contigs.fasta", output=OUTPUT, sample=set(EXPS["Sample"])),
        expand("{output}/Annotation/{sample}/fgs.ffn", output=OUTPUT, sample=set(EXPS["Sample"]))
    output:
        expand("{output}/Quantification/{name}.readcounts", output=OUTPUT, name=set(not_mp_exps['Name'])),
        expand("{output}/Quantification/{sample}_mg.readcounts", output=OUTPUT, sample=set(mg_exps['Sample'])),
        expand("{output}/Quantification/{sample}_mg_norm.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])),
        expand("{output}/Quantification/{sample}_mt.readcounts", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Quantification/{sample}_mt_norm.tsv", output=OUTPUT, sample=set(mg_exps['Sample']))
    threads:
        config["threads"]
    params:
        output = OUTPUT,
        exps = f"{OUTPUT}/exps.tsv"
    conda:
        "../envs/quantification.yaml"
    script:
        "../scripts/quantification.py"
