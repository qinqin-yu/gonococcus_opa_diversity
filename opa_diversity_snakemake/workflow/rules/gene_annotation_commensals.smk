rule prokka_commensals:
    input:
        fasta="input_data/complete_genome_assemblies/commensals/{sample}.fa"
    output:
        gff="results/commensals/annotations/{sample}.gff"
    params:
        name="{sample}",
        tmpdir=SCRATCH+"/commensals/annotations/{sample}"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
        time=lambda wildcards, attempt: attempt * 10,
        cpus=8
    threads: 8
    conda:
        "../envs/prokka.yml"
    shell:
        """
        mkdir -p results/commensals/annotations/
        prokka --force --centre X --compliant --outdir {params.tmpdir} --locustag {params.name} --prefix {params.name} --genus Neisseria --strain {params.name} --cpus 8 {input.fasta}
        cp {params.tmpdir}/*.gff results/commensals/annotations
        rm -r {params.tmpdir}
        """

rule parse_prokka_commensals:
    input:
        gff="results/commensals/annotations/{sample}.gff"
    output:
        csv="results/commensals/annotations/{sample}.csv"
    conda:
        "../envs/bcbio-gff.yml"
    shell:
        """
        workflow/scripts/prokka_gff2csv.py {input.gff} {output.csv}
        """