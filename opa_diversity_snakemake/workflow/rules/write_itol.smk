rule get_opa_num:
    input:
        "results/opa_metadata_locus.csv"
    output:
        "results/summary_statistics/num_opa_genes_by_strain.csv"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            workflow/scripts/get_opa_num.py {input} {output}
        """

rule write_itol_num_opa_genes:
    input:
        "results/summary_statistics/num_opa_genes_by_strain.csv"
    output:
        "results/itol/itol_num_opa_genes.txt"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            mkdir -p results/itol
            workflow/scripts/write_itol_num_opa_genes.py {input} {output}
        """

rule write_itol_strain_metadata:
    input:
        opa_metadata="results/opa_metadata_locus.csv",
        genome_metadata="input_data/lab_strains_hybrid_genomes_metadata.csv"
    output:
        anatomic_site="results/itol/itol_anatomic_site.txt",
        gender="results/itol/itol_gender.txt",
        sexual_behavior="results/itol/itol_sexual_behavior.txt"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            mkdir -p results/itol
            workflow/scripts/write_itol_strain_metadata.py {input.opa_metadata} {input.genome_metadata} {output.anatomic_site} {output.gender} {output.sexual_behavior}
        """