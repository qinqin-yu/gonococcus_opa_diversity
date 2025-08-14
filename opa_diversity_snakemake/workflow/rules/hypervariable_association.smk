rule pseudogenome_alignment_fastbaps:
    input:
        "input_data/complete_genome_pseudogenomes.fasta.all"
    output:
        "results/fastbaps/complete_genome_pseudogenomes_renamed.fasta.all"
    conda:
        "../envs/biopython_plotting.yml"
    shell:
        """
            mkdir -p results/fastbaps/
            workflow/scripts/get_pseudogenomes_for_fastbaps.py {input} {output}
        """

rule fastbaps:
    input:
        pseudogenome_alignment="results/fastbaps/complete_genome_pseudogenomes_renamed.fasta.all",
        tree="results/gubbins/complete_genome_pseudogenomes.final_tree.tre"
    output:
        clusters="results/fastbaps/fastbaps_clusters.csv"
    conda:
        "../envs/fastbaps_argparse.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        time=lambda wildcards, attempt: attempt * 60
    shell:
        """
            mkdir -p fastbaps/
            workflow/scripts/fastbaps_whole_genome.R --pseudogenome_alignment {input.pseudogenome_alignment} --tree {input.tree} --clusters {output.clusters}
        """

rule hv1_hv2_association:
    input:
        fastbaps_clusters='results/fastbaps/fastbaps_clusters.csv',
        mcl_clusters='results/mcl/mcl_clusters.csv'
    output:
        "results/fastbaps/hv1_hv2_association_results.txt"
    conda:
        "../envs/munkres.yml"
    shell:
        """
            workflow/scripts/hv1_hv2_association_subsample_baps.py {input.fastbaps_clusters} {input.mcl_clusters} {output}
        """

