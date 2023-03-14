rule preprocess:
    input:
        lambda wildcards: EXPS.loc[EXPS['Name'] == wildcards.name, 'Files'].iloc[0].split(',')
    output:
        expand("{output}/Preprocess/Trimmomatic/quality_trimmed_{{name}}{fr}.fq", output=OUTPUT,
            fr=(['_forward_paired', '_reverse_paired'] if EXPS["Files"].str.contains(',').tolist() else ''))
    threads:
        config["threads"]
    params:
        reads = lambda wildcards: EXPS.loc[EXPS['Name'] == wildcards.name, 'Files'].iloc[0],
        resources_directory = config["resources_directory"],
        data_type = lambda wildcards: EXPS.loc[EXPS['Name'] == wildcards.name, 'Data type'].iloc[0],
        minlen = config["minimum_read_length"],
        avgqual = config["minimum_read_average_quality"]
    conda:
        "../envs/preprocess.yaml"
    shell:
        "python ../scripts/preprocess.py -i {params.reads} -t {threads} -o {OUTPUT}/Preprocess "
        "-d {params.data_type} -rd {params.resources_directory} -n {wildcards.name} --minlen {params.minlen} "
        "--avgqual {params.avgqual}"
