rule keggcharter:
    input:
        f"{OUTPUT}/MOSCA_Entry_Report.tsv"
    output:
        f"{OUTPUT}/KEGG_maps/KEGGCharter_results.tsv"
    threads:
        1
    params:
        outdir=f"{OUTPUT}/KEGG_maps",
        resources_directory=config["resources_directory"],
        metabolic_maps=f" -mm {','.join(config['keggcharter_maps']) if len(config['keggcharter_maps']) > 0 else ''}",
        quant_part=f"-qcol {','.join(not_mg_exps['Name'].tolist())}" if len(not_mg_exps) > 0 else "-iq",
        taxa_level=config["keggcharter_taxa_level"],
        number_of_taxa=config["keggcharter_number_of_taxa"]
    conda:
        "../envs/keggcharter.yaml"
    shell:
        "keggcharter -f {input} -o {params.outdir} {params.quant_part} -tc 'Taxonomic lineage ({params.taxa_level})' "
        "-not {params.number_of_taxa} -keggc KEGG -ecc 'EC number' -cogc 'COG ID' "     # TODO - missing KO column
        "-rd {params.resources_directory}{params.metabolic_maps}"
