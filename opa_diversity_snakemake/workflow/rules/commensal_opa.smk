rule opa_id_commensals:
    # Note that these complete genomes are not rotated
    input:
        "input_data/complete_genome_assemblies/commensals/{sample}.fa"
    output:
        "results/commensals/opa_sequences/{sample}.fa",
        "results/commensals/opa_sequences_no_repeats/{sample}.fa",
        "results/commensals/opa_locations/{sample}.csv"
    conda:
        "../envs/biopython_regex.yml"
    params:
        outdir="results/commensals/"
    shell:
        """
            mkdir -p results/commensals/opa_sequences/
            mkdir -p results/commensals/opa_sequences_no_repeats/
            mkdir -p results/commensals/opa_locations/
            workflow/scripts/identify_opa_genes.py {input} {params.outdir}
        """

rule postprocess_commensal_opa:
    input:
        commensal_metadata="input_data/complete_genome_assemblies/commensals/commensals_metadata.csv"
    output:
        opa_commensals_metadata="results/commensals/opa_commensals_metadata.csv"
    params:
        opa_locations_folder="results/commensals/opa_locations/"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            workflow/scripts/postprocess_commensal_opa.py {input.commensal_metadata} {params.opa_locations_folder} {output.opa_commensals_metadata}
        """

rule mash_neisseria:
    input:
        commensals=expand("results/commensals/opa_sequences_no_repeats/{sample}.fa", sample = SAMPLES_COMMENSALS),
        gc="results/opa_sequences_no_repeats/opa_sequences_no_repeats.fa.all"
    output:
        neisseria="results/commensals/mash/opa_sequences_no_repeats_neisseria.fa",
        sketch="results/commensals/mash/opa_sequences_no_repeats_neisseria.fa.msh",
        dist="results/commensals/mash/neisseria_mash_distance.tab"
    conda:
        "../envs/mash.yml"
    shell:
        """
            mkdir -p mash/

            cat {input.commensals}>results/commensals/opa_sequences_no_repeats/opa_sequences_no_repeats.fa.all
            cat results/commensals/opa_sequences_no_repeats/opa_sequences_no_repeats.fa.all {input.gc} > {output.neisseria}
            
            # Sketching
            mash sketch -i -k 9 {output.neisseria}

            # Calculate distance matrix
            mash dist {output.sketch} {output.sketch} > {output.dist}
        """

rule mash_to_phylip_neisseria:
    input:
        mash="results/commensals/mash/neisseria_mash_distance.tab"
    output:
        phylip="results/commensals/mash/neisseria_mash_distance.phy"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            workflow/scripts/convert_mash_distance_matrix_to_phylip.py {input.mash} {output.phylip}
        """

rule rapidnj_neisseria:
    input:
        phylip="results/commensals/mash/neisseria_mash_distance.phy"
    output:
        tree="results/commensals/rapidnj/neisseria_opa.nwk"
    conda:
        "../envs/rapidnj.yml"
    shell:
        """
            mkdir -p commensals/rapidnj/
            rapidnj {input.phylip} -i pd -x {output.tree}

            # Remove apostrophes from opa names in tree
            sed -i "s/[']//g" {output.tree}
        """

rule itol_neisseria:
    input:
        neisseria="results/commensals/opa_commensals_metadata.csv",
        gc="results/opa_metadata_locus.csv"
    output:
        "results/commensals/itol/itol_neisseria_species.txt"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            mkdir -p results/commensals/itol/
            workflow/scripts/write_itol_neisseria.py {input.neisseria} {input.gc} {output}
        """

rule mafft_neisseria_opa:
    input:
        "results/commensals/mash/opa_sequences_no_repeats_neisseria.fa"
    output:
        "results/commensals/mafft/opa_sequences_no_repeats_neisseria.fa.aln"
    conda:
        "../envs/mafft.yml"
    shell:
        """
            mkdir -p results/commensals/mafft
            mafft --auto {input} > {output}
        """