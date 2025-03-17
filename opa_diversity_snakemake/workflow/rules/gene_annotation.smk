rule prokka:
    input:
        fasta="results/assemblies_shifted/{sample}.fa"
    output:
        gff="results/annotations/{sample}.gff"
    params:
        name="{sample}",
        tmpdir=SCRATCH+"/annotations/{sample}"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
        time=lambda wildcards, attempt: attempt * 10,
        cpus=8
    threads: 8
    conda:
        "../envs/prokka.yml"
    shell:
        """
        mkdir -p results/annotations/
        prokka --force --centre X --compliant --outdir {params.tmpdir} --locustag {params.name} --prefix {params.name} --genus Neisseria --species gonorrhoeae --strain {params.name} --cpus 8 {input.fasta}
        cp {params.tmpdir}/*.gff results/annotations
        rm -r {params.tmpdir}
        """

rule parse_prokka:
    input:
        gff="results/annotations/{sample}.gff"
    output:
        csv="results/annotations/{sample}.csv"
    conda:
        "../envs/bcbio-gff.yml"
    shell:
        """
        workflow/scripts/prokka_gff2csv.py {input.gff} {output.csv}
        """