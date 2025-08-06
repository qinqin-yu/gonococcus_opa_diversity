rule simulate_short_reads_public:
    input:
        "input_data/complete_genome_assemblies/public/{sample}.fa"
    output:
        expand(SCRATCH+"/simulated_short_reads/{{sample}}_{paired}.fastq.gz", paired=PAIRED)
    params:
        tmpdir=SCRATCH+"/simulated_short_reads"
    conda:
        "../envs/art.yml"
    shell:
        """
            mkdir -p {params.tmpdir}
            art_illumina -p -ss MSv3 -i {input} -l 250 -f 80 -m 600 -s 100 -o {params.tmpdir}/{wildcards.sample}_ -na
            mv {params.tmpdir}/{wildcards.sample}_1.fq {params.tmpdir}/{wildcards.sample}_1.fastq
            mv {params.tmpdir}/{wildcards.sample}_2.fq {params.tmpdir}/{wildcards.sample}_2.fastq
            gzip {params.tmpdir}/{wildcards.sample}_1.fastq
            gzip {params.tmpdir}/{wildcards.sample}_2.fastq
        """

rule simulate_short_reads_public_new:
    input:
        "input_data/complete_genome_assemblies/public_20230215_20250731/{sample}.fa"
    output:
        expand(SCRATCH+"/simulated_short_reads/{{sample}}_{paired}.fastq.gz", paired=PAIRED)
    params:
        tmpdir=SCRATCH+"/simulated_short_reads"
    conda:
        "../envs/art.yml"
    shell:
        """
            mkdir -p {params.tmpdir}
            art_illumina -p -ss MSv3 -i {input} -l 250 -f 80 -m 600 -s 100 -o {params.tmpdir}/{wildcards.sample}_ -na
            mv {params.tmpdir}/{wildcards.sample}_1.fq {params.tmpdir}/{wildcards.sample}_1.fastq
            mv {params.tmpdir}/{wildcards.sample}_2.fq {params.tmpdir}/{wildcards.sample}_2.fastq
            gzip {params.tmpdir}/{wildcards.sample}_1.fastq
            gzip {params.tmpdir}/{wildcards.sample}_2.fastq
        """

rule simulate_short_reads_grad_lab:
    input:
        "input_data/complete_genome_assemblies/grad_lab/{sample}.fa"
    output:
        expand(SCRATCH+"/simulated_short_reads/{{sample}}_{paired}.fastq.gz", paired=PAIRED)
    params:
        tmpdir=SCRATCH+"/simulated_short_reads"
    conda:
        "../envs/art.yml"
    shell:
        """
            mkdir -p {params.tmpdir}
            art_illumina -p -ss MSv3 -i {input} -l 250 -f 80 -m 600 -s 100 -o {params.tmpdir}/{wildcards.sample}_ -na
            mv {params.tmpdir}/{wildcards.sample}_1.fq {params.tmpdir}/{wildcards.sample}_1.fastq
            mv {params.tmpdir}/{wildcards.sample}_2.fq {params.tmpdir}/{wildcards.sample}_2.fastq
            gzip {params.tmpdir}/{wildcards.sample}_1.fastq
            gzip {params.tmpdir}/{wildcards.sample}_2.fastq
        """

rule gubbins:
    input:
        expand("input_data/complete_genome_pseudogenomes/{sample}_pseudogenome.fasta", sample=SAMPLES)
    output:
        protected("results/gubbins/complete_genome_pseudogenomes.filtered_polymorphic_sites.fasta"),
        protected("results/gubbins/complete_genome_pseudogenomes.final_tree.tre")
    params:
        tmpdir=SCRATCH+"/gubbins"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 12000,
        time=lambda wildcards, attempt: attempt * 360,
        cpus=12
    conda:
        "../envs/gubbins.yml"
    shell:
        """
            mkdir -p {params.tmpdir} #results/gubbins
            cat {input} > {params.tmpdir}/complete_genome_pseudogenomes.fasta #results/gubbins/complete_genome_pseudogenomes.fasta
            run_gubbins.py --first-tree-builder rapidnj --tree-builder raxmlng --first-model JC --model GTR --threads 12 --iterations 20 --prefix {params.tmpdir}/complete_genome_pseudogenomes {params.tmpdir}/complete_genome_pseudogenomes.fasta #results/gubbins/complete_genome_pseudogenomes.fasta

            mkdir -p results/gubbins
            cp {params.tmpdir}/complete_genome_pseudogenomes.filtered_polymorphic_sites.fasta results/gubbins
            cp {params.tmpdir}/complete_genome_pseudogenomes.final_tree.tre results/gubbins
        """

rule remove_duplicate_sequences:
    input:
        expand("input_data/complete_genome_pseudogenomes/{sample}_pseudogenome.fasta", sample=SAMPLES),
        expand("input_data/representative_genome_pseudogenomes/{sample}_pseudogenome.fasta", sample=SAMPLES_REPRESENTATIVE)
    output:
        SCRATCH+"/gubbins_representative_and_complete/representative_and_complete_genome_pseudogenomes_dedup.fasta"
    params:
        tmpdir=SCRATCH+"/gubbins_representative_and_complete"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
            mkdir -p {params.tmpdir}
            cat {input} > {params.tmpdir}/representative_and_complete_genome_pseudogenomes.fasta
            seqkit rmdup {params.tmpdir}/representative_and_complete_genome_pseudogenomes.fasta>{params.tmpdir}/representative_and_complete_genome_pseudogenomes_dedup.fasta
        """

rule gubbins_representative_and_complete_genomes:
    input:
        SCRATCH+"/gubbins_representative_and_complete/representative_and_complete_genome_pseudogenomes_dedup.fasta"
    output:
        protected("results/gubbins/representative_and_complete_genome_pseudogenomes_dedup.filtered_polymorphic_sites.fasta"),
        protected("results/gubbins/representative_and_complete_genome_pseudogenomes_dedup.final_tree.tre")
    params:
        tmpdir=SCRATCH+"/gubbins_representative_and_complete"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 12000,
        time=10080,
        cpus=12,
        partition="intermediate"
    conda:
        "../envs/gubbins.yml"
    shell:
        """
            run_gubbins.py --first-tree-builder rapidnj --tree-builder raxmlng --first-model JC --model GTR --threads 12 --prefix {params.tmpdir}/representative_and_complete_genome_pseudogenomes {params.tmpdir}/representative_and_complete_genome_pseudogenomes_dedup.fasta

            mkdir -p results/gubbins
            cp {params.tmpdir}/representative_and_complete_genome_pseudogenomes_dedup.filtered_polymorphic_sites.fasta results/gubbins
            cp {params.tmpdir}/representative_and_complete_genome_pseudogenomes_dedup.final_tree.tre results/gubbins
        """