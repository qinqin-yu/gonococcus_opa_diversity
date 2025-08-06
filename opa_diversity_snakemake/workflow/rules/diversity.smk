rule translate_opa_seqs:
    input:
        "results/opa_sequences_no_repeats/{sample}.fa"
    output:
        "results/opa_sequences_aa/{sample}.fa"
    conda:
        "../envs/biopython_regex.yml"
    shell:
        """
            mkdir -p "results/opa_sequences_aa/"
            workflow/scripts/get_opa_orfs_aa.py {input} {output}
        """

rule filter_opa_aa_seqs:
    input:
        expand("results/opa_sequences_aa/{sample}.fa", sample=SAMPLES)
    output:
        fasta="results/opa_sequences_aa/opa_sequences_aa.fa.all",
        fasta_clean="results/opa_sequences_aa/opa_sequences_aa.fa.all.clean",
        excluded_sequences="results/opa_sequences_aa/opa_sequences_excluded_in_alignment.csv"
    conda:
        "../envs/biopython_regex.yml"
    shell:
        """
            cat {input} > {output.fasta}
            workflow/scripts/get_opa_orfs_aa_for_alignment.py {output.fasta} {output.fasta_clean} {output.excluded_sequences}
        """

rule mafft:
    input:
        "results/opa_sequences_aa/opa_sequences_aa.fa.all.clean"
    output:
        "results/alignments/opa_sequences_aa.fa.aln"
    conda:
        "../envs/mafft.yml"
    resources:
        mem_mb=4000
    shell:
        """
            mkdir -p "results/alignments/"
            mafft --auto {input} > {output}
        """

rule calculate_pairwise_distances:
    input:
        "results/alignments/opa_sequences_aa.fa.aln"
    output:
        pairwise_distances="results/pairwise_distance/within_strain_pairwise_distance_aa.csv",
        png="figures/diversity/pairwise_distance_aa_within_between_strains_hist.png",
        pdf="figures/diversity/pairwise_distance_aa_within_between_strains_hist.pdf"
    conda:
        "../envs/biopython_plotting.yml"
    shell:
        """
            mkdir -p results/pairwise_distance
            mkdir -p figures/diversity
            workflow/scripts/calculate_pairwise_distances_aa.py {input} {output.pairwise_distances} {output.png} {output.pdf}
        """

rule similar_opa:
    input:
        'results/pairwise_distance/within_strain_pairwise_distance_aa.csv'
    output:
        similar_opas='results/pairwise_distance/within_strain_similar_opa_aa.csv',
        summary_png='figures/diversity/within_strain_similar_opa_aa.png',
        summary_pdf='figures/diversity/within_strain_similar_opa_aa.pdf',
        vary_thresh_png='figures/diversity/within_strain_similar_opa_aa_num_strains.png',
        vary_thresh_pdf='figures/diversity/within_strain_similar_opa_aa_num_strains.pdf'
    conda:
        "../envs/networkx_plotting.yml"
    shell:
        """
            mkdir -p figures/diversity
            workflow/scripts/within_strain_similar_opa_aa.py {input} {output.similar_opas} {output.summary_png} {output.summary_pdf} {output.vary_thresh_png} {output.vary_thresh_pdf}
        """

rule pseudogenome_distances:
    input:
        expand("input_data/complete_genome_pseudogenomes/{sample}_pseudogenome.fasta", sample = SAMPLES)
    output:
        "results/pairwise_distance/pseudogenome_distances.csv"
    conda:
        "../envs/pairsnp.yml"
    shell:
        """
        cat {input}>input_data/complete_genome_pseudogenomes.fasta.all
        pairsnp -c input_data/complete_genome_pseudogenomes.fasta.all > {output}
        """

rule opa_repertoire_distance:
    input:
        alignment="results/alignments/opa_sequences_aa.fa.aln",
        phylogeny="results/itol/complete_genome_pseudogenomes.final_tree.itol_order.tre",
        pairsnp_output="results/pairwise_distance/pseudogenome_distances.csv"
    output:
        distances = 'results/pairwise_distance/between_strain_opa_distance_aa.csv',
        summary = 'results/pairwise_distance/between_strain_opa_and_genome_distance_aa.csv'
    conda:
        "../envs/biopython_plotting.yml"
    shell:
        """
            workflow/scripts/calculate_opa_repertoire_distance.py {input.alignment} {input.phylogeny} {input.pairsnp_output} {output.distances} {output.summary}
        """