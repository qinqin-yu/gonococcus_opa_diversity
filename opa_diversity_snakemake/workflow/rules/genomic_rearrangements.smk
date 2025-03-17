rule progressivemauve:
    input:
        assembly="results/assemblies_shifted/{sample}.fa",
        reference="results/assemblies_shifted/FA1090.fa"
    output:
        xmfa="results/progressivemauve/pairwise_outputs/{sample}.xmfa"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 12000,
        time=lambda wildcards, attempt: attempt * 60
    conda:
        "../envs/mauve.yml"
    shell:
        """
            mkdir -p results/progressivemauve/pairwise_outputs
            progressiveMauve --output={output.xmfa} {input.reference} {input.assembly}
        """

rule parse_progressivemauve:
    input:
        xmfa=expand("results/progressivemauve/pairwise_outputs/{sample}.xmfa", sample=SAMPLES),
        opa_location=expand("results/opa_locations/{sample}.csv", sample=SAMPLES)
    output:
        opa_metadata=temp("results/opa_metadata.csv"),
        lcbs="results/progressivemauve/pairwise_lcbs.csv"
    params:
        xmfa_dir="results/progressivemauve/pairwise_outputs/",
        opa_locations_dir="results/opa_locations/"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            workflow/scripts/parse_progressivemauve.py {params.xmfa_dir} {params.opa_locations_dir} {output.lcbs} {output.opa_metadata}
        """

rule assign_opa_loci:
    input:
        "results/opa_metadata.csv"
    output:
        metadata="results/opa_metadata_locus.csv",
        locus_assignment="results/progressivemauve/locus_assignment_ranges.csv",
        fig_filename="figures/opa_loci_definitions.png"

    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            workflow/scripts/assign_opa_loci.py {input} {output.metadata} {output.locus_assignment} {output.fig_filename}
        """