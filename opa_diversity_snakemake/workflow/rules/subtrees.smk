rule mafft_opa:
    input:
        unaligned_opa = "results/subtrees/opa_sequences_collated/{sample}.fa",
    output: 
        aligned_opa = "results/subtrees/opa_sequences_aligned/{sample}.aln"
    conda:
        "../envs/mafft.yml"
    shell:
        """
        mkdir -p results/subtrees/opa_sequences_aligned
        mafft --auto {input.unaligned_opa}>{output.aligned_opa}
        """

# Maybe consider adding rules later that determine which reference genomes to use for subtrees (currently in ipynb notebook)

def get_reference(sample):
    # Get the subtree reference to map to
    df = pd.read_csv("results/subtrees/pseudogenome_references/pseudogenome_reference.csv")
    reference = df[df['isolate']==sample]['reference'].values[0]
    reference_path = "results/subtrees/pseudogenome_references/"+reference+".fa"
    return reference_path

def get_indexes(sample):
    # Get the subtree reference to map to
    df = pd.read_csv("results/subtrees/pseudogenome_references/pseudogenome_reference.csv")
    reference = df[df['isolate']==sample]['reference'].values[0]
    indexes = ["results/subtrees/pseudogenome_references/"+reference+".fa.amb",
                "results/subtrees/pseudogenome_references/"+reference+".fa.ann",
                "results/subtrees/pseudogenome_references/"+reference+".fa.bwt",
                "results/subtrees/pseudogenome_references/"+reference+".fa.pac",
                "results/subtrees/pseudogenome_references/"+reference+".fa.sa"]
    return indexes

rule subtree_reference_bwa_index:
    input:
        "results/subtrees/pseudogenome_references/{sample}.fa"
    output:
        multiext("results/subtrees/pseudogenome_references/{sample}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    wrapper:
        "v5.10.0/bio/bwa/index"
    
rule subtree_map_sort:
    input:
        reads=expand(SCRATCH+"/simulated_short_reads/{{sample}}_{paired}.fastq.gz", paired=PAIRED),
        idx=lambda wildcards: get_indexes(wildcards.sample)
    output:
        temp(SCRATCH+"/subtree_mapping/{sample}.bam")
    params:
        index=lambda wildcards: get_reference(wildcards.sample),
        sorting="samtools",
        sort_order="coordinate"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
        time=lambda wildcards, attempt: attempt * 35
    wrapper:
        "v1.28.0/bio/bwa/mem"

rule subtree_duplicates:
    input:
        bams=SCRATCH+"/subtree_mapping/{sample}.bam"
    output:
        bam=SCRATCH+"/subtree_mapping/{sample}.marked.bam",
        metrics=temp(SCRATCH+"/subtree_mapping/{sample}.metrics")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000,
        time=lambda wildcards, attempt: attempt * 10
    wrapper:
        "v1.28.0/bio/picard/markduplicates"

rule subtree_index:
    input:
        SCRATCH+"/subtree_mapping/{sample}.marked.bam"
    output:
        SCRATCH+"/subtree_mapping/{sample}.marked.bam.bai"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 100,
        time=lambda wildcards, attempt: attempt * 1
    wrapper:
        "v1.28.0/bio/samtools/index"

rule subtree_variants:
    input:
        bam=SCRATCH+"/subtree_mapping/{sample}.marked.bam",
        index=SCRATCH+"/subtree_mapping/{sample}.marked.bam.bai",
        reference=lambda wildcards: get_reference(wildcards.sample)
    output:
        SCRATCH+"/subtree_variants/{sample}_pilon.vcf.gz",
        temp(SCRATCH + "/subtree_variants/{sample}_pilon.fasta")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
        time=lambda wildcards, attempt: attempt * 60
    params:
        tmpdir=SCRATCH+'/subtree_variants'
    conda:
        "../envs/tabix.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        java -Xmx16g -jar /n/holylfs05/LABS/grad_lab/Lab/software/pilon-1.24.jar --genome {input.reference} --bam {input.bam} --output {params.tmpdir}/{wildcards.sample}_pilon --variant --vcf --mindepth 10 --minmq 20
        bgzip {params.tmpdir}/{wildcards.sample}_pilon.vcf
        """

rule subtree_pseudogenome:
    input:
        SCRATCH+"/subtree_variants/{sample}_pilon.vcf.gz"
    output:
        "results/subtrees/pseudogenomes/{sample}_subtree_pseudogenome.fasta"
    params:
        tmpdir=SCRATCH+"/subtree_pseudogenomes"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 400, # Change back to 100 later
        time=lambda wildcards, attempt: attempt * 1
    conda:
        "../envs/biopython_plotting.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        workflow/scripts/pilonVCFtoFasta_AF0.9_QY.py {input} {params.tmpdir}
        cp {params.tmpdir}/{wildcards.sample}_pseudogenome.fasta {output}
        """

def get_isolates_in_subtree(subtree):
    # Get the isolates in the subtree
    isolates = pd.read_csv("results/subtrees/isolates/subtree_" + str(subtree) + "_isolates.txt", header = None)[0].values
    pseudogenome_paths = []
    for isolate in isolates:
        pseudogenome_paths.append("results/subtrees/pseudogenomes/" + str(isolate) + "_subtree_pseudogenome.fasta")
    print(pseudogenome_paths)
    return pseudogenome_paths

rule subtree_gubbins:
    input:
        lambda wildcards: get_isolates_in_subtree(wildcards.subtree)
    output:
        protected("results/subtrees/gubbins/subtree_{subtree}.filtered_polymorphic_sites.fasta"),
        protected("results/subtrees/gubbins/subtree_{subtree}.final_tree.tre")
    params:
        tmpdir=SCRATCH+"/subtree_gubbins"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 12000,
        time=lambda wildcards, attempt: attempt * 360,
        cpus=12
    conda:
        "../envs/gubbins.yml"
    shell:
        """
            mkdir -p {params.tmpdir}
            cat {input} > {params.tmpdir}/subtree_{wildcards.subtree}.fasta
            run_gubbins.py --first-tree-builder rapidnj --tree-builder raxmlng --first-model JC --model GTR --threads 12 --iterations 20 --prefix {params.tmpdir}/subtree_{wildcards.subtree} {params.tmpdir}/subtree_{wildcards.subtree}.fasta

            mkdir -p results/subtrees/gubbins
            cp {params.tmpdir}/subtree_{wildcards.subtree}.filtered_polymorphic_sites.fasta results/subtrees/gubbins
            cp {params.tmpdir}/subtree_{wildcards.subtree}.final_tree.tre results/subtrees/gubbins
        """

rule subtree_get_inputs_for_beast:
    input:
        metadata='input_data/lab_strains_hybrid_genomes_metadata.csv',
        alignment='results/subtrees/gubbins/subtree_{subtree}.filtered_polymorphic_sites.fasta'
    output:
        alignment_with_dates='results/subtrees/beast/subtree_{subtree}_recomb_masked_with_dates.fasta',
        dates='results/subtrees/beast/subtree_{subtree}_dates.dat'
    conda:
        '../envs/biopython_plotting.yml'
    shell:
        """
            mkdir -p results/subtrees/beast
            workflow/scripts/get_inputs_for_beast.py {input.metadata} {input.alignment} {output.alignment_with_dates} {output.dates}
        """

# rule subtree_beast:
#     input:
#         'results/subtrees/beast/subtree_{subtree}_recomb_masked_with_dates.xml'
#     output:
#         tree="results/subtrees/beast/subtree_{subtree}_recomb_masked_with_dates-subtree_{subtree}_recomb_masked_with_dates.trees",
#         log="results/subtrees/beast/subtree_{subtree}_recomb_masked_with_dates.log"
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 12000,
#         time=lambda wildcards, attempt: attempt * 360,
#         cpus=12
#     shell:
#         """
#             /n/holylfs05/LABS/grad_lab/Users/qinqinyu/software/beast/bin/beast -threads 12 -working
#         """
