rule keggcharter:
    input:
        f"{OUTPUT}/MOSCA_Entry_Report.xlsx"
    output:
        f"{OUTPUT}/KEGG_maps/KEGGCharter_results.tsv"
    threads:
        1
    params:
        outdir=f"{OUTPUT}/KEGG_maps",
        mg_cols= ','.join(mg_exps['Name'].tolist()),
        resources_directory=config["resources_directory"],
        metabolic_maps=f" -mm {','.join(config['keggcharter_maps']) if len(config['keggcharter_maps']) > 0 else ''}",
        exp_cols=f" -tcol {','.join(mt_exps['Name'].tolist())}" if len(mt_exps) > 0 else '',
        taxa_level=config["keggcharter_taxa_level"],
        number_of_taxa=config["keggcharter_number_of_taxa"]
    conda:
        "../envs/keggcharter.yaml"
    shell:
        "keggcharter.py -f {input} -o {params.outdir} -gcol {params.mg_cols}{params.exp_cols} "
        "-tc 'Taxonomic lineage ({params.taxa_level})' -not {params.number_of_taxa} -keggc KEGG "
        "-rd {params.resources_directory}{params.metabolic_maps}"
