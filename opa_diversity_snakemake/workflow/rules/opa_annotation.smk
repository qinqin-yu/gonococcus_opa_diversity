rule rotate_origin_grad_lab:
    input:
        "input_data/complete_genome_assemblies/grad_lab/{sample}.fa",
    output: 
        "results/assemblies_shifted/{sample}.fa"
    shell:
        """
        mkdir -p results/assemblies_shifted/
        # Match to first 90 bases of dnaA in FA1090 reference sequence (NC_002946.2) allowing for up to 5 mismatches (mening has a very similar sequence)
        workflow/scripts/rotate/rotate -s ATGACATTAGCAGAGTTTTGGCCGCTGTGCCTCCGCCGTCTTCACGATATGTTGCCTCACGGGCAGTTTGCGCAATGGATTGCGCCCCTT -m 5 {input} > {output}
        """

rule rotate_origin_public:
    input:
        "input_data/complete_genome_assemblies/public/{sample}.fa",
    output: 
        "results/assemblies_shifted/{sample}.fa"
    shell:
        """
        mkdir -p results/assemblies_shifted/
        # Match to first 90 bases of dnaA in FA1090 reference sequence (NC_002946.2) allowing for up to 5 mismatches (mening has a very similar sequence)
        workflow/scripts/rotate/rotate -s ATGACATTAGCAGAGTTTTGGCCGCTGTGCCTCCGCCGTCTTCACGATATGTTGCCTCACGGGCAGTTTGCGCAATGGATTGCGCCCCTT -m 5 {input} > {output}
        """

rule rotate_origin_public_new:
    input:
        "input_data/complete_genome_assemblies/public_20230215_20250731/{sample}.fa",
    output: 
        "results/assemblies_shifted/{sample}.fa"
    shell:
        """
        mkdir -p results/assemblies_shifted/
        # Match to first 90 bases of dnaA in FA1090 reference sequence (NC_002946.2) allowing for up to 5 mismatches (mening has a very similar sequence)
        workflow/scripts/rotate/rotate -s ATGACATTAGCAGAGTTTTGGCCGCTGTGCCTCCGCCGTCTTCACGATATGTTGCCTCACGGGCAGTTTGCGCAATGGATTGCGCCCCTT -m 5 {input} > {output}
        """

rule rotate_origin_public_new_2022NG_0032:
    # Note that this is because the starting sequence has a 1 bp deletion so that only the sequence ATGACATTAGCAGAGTTTT will produce a match
    input:
        "input_data/complete_genome_assemblies/public_20230215_20250731/2022NG-0032.fa",
    output: 
        "results/assemblies_shifted/2022NG-0032.fa"
    shell:
        """
        mkdir -p results/assemblies_shifted/
        # Match to first 19 bases of dnaA in FA1090 reference sequence (NC_002946.2) allowing for up to 0 mismatches
        workflow/scripts/rotate/rotate -s ATGACATTAGCAGAGTTTT {input} > {output}
        """

rule opa_id:
    input:
        "results/assemblies_shifted/{sample}.fa"
    output:
        "results/opa_sequences/{sample}.fa",
        "results/opa_sequences_no_repeats/{sample}.fa",
        "results/opa_locations/{sample}.csv"
    conda:
        "../envs/biopython_regex.yml"
    params:
        outdir="results/"
    shell:
        """
            mkdir -p results/opa_sequences/
            mkdir -p results/opa_sequences_no_repeats/
            mkdir -p results/opa_locations/
            workflow/scripts/identify_opa_genes.py {input} {params.outdir}
        """