rule preprocess:
    input:
        lambda wildcards: EXPS.loc[EXPS['Name'] == wildcards.name, 'Files'].iloc[0].split(',')
    output:
        expand("{output}/Preprocess/Trimmomatic/quality_trimmed_{{name}}{fr}.fq", output=OUTPUT,
            fr=(['_forward_paired', '_reverse_paired'] if EXPS["Files"].str.contains(',').tolist() else ''))
    threads:
        config["threads"]
    params:
        output = f'{OUTPUT}/Preprocess',
        name = lambda wildcards: wildcards.name,
        reads = lambda wildcards: EXPS.loc[EXPS['Name'] == wildcards.name, 'Files'].iloc[0],
        resources_directory = config["resources_directory"],
        data_type = lambda wildcards: EXPS.loc[EXPS['Name'] == wildcards.name, 'Data type'].iloc[0],
        rrna_db = config["sortmerna_database"],
        mg_minlen = config["minimum_mg_read_length"],
        mt_minlen= config["minimum_mt_read_length"],
        avgqual = config["minimum_read_average_quality"]
    conda:
        "../envs/preprocess.yaml"
    script:
        "../scripts/preprocess.py"
