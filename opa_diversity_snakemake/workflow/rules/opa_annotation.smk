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