rule mafft_nt:
    input:
        expand("results/opa_sequences_no_repeats/{sample}.fa", sample=SAMPLES)
    output:
        fasta="results/opa_sequences_no_repeats/opa_sequences_no_repeats.fa.all",
        aln="results/alignments/opa_sequences_no_repeats.fa.aln"
    params:
        ref="resources/FA1090_opa_1_reference.fa"
    conda:
        "../envs/mafft.yml"
    resources:
        mem_mb=4000
    shell:
        """
            mkdir -p "results/alignments/"
            cat {input} {params.ref} > {output.fasta}
            # Note that even though muscle may be more accurate, it is extremely slow on this number of sequences. Thus, I have decided to use mafft with a manually adjusted gap opening penalty (--op) and gap extension penalty (--ep). 
            mafft --auto --op 4 --ep 1 {output.fasta} > {output.aln}
        """

rule split_alignment:
    input:
        "results/alignments/opa_sequences_no_repeats.fa.aln"
    output:
        sv="results/alignments/opa_sequences_no_repeats_sv.fa.aln",
        hv1="results/alignments/opa_sequences_no_repeats_hv1.fa.aln",
        hv2="results/alignments/opa_sequences_no_repeats_hv2.fa.aln",
        sv_ungapped="results/mash/opa_sequences_no_repeats_sv_ungapped.fa",
        hv1_ungapped="results/mash/opa_sequences_no_repeats_hv1_ungapped.fa",
        hv2_ungapped="results/mash/opa_sequences_no_repeats_hv2_ungapped.fa"
    conda:
        "../envs/biopython_plotting.yml"
    shell:
        """
            mkdir -p results/mash
            workflow/scripts/split_alignment.py {input} {output.sv} {output.hv1} {output.hv2} {output.sv_ungapped} {output.hv1_ungapped} {output.hv2_ungapped}
        """

rule mash:
    input:
        sv_ungapped="results/mash/opa_sequences_no_repeats_sv_ungapped.fa",
        hv1_ungapped="results/mash/opa_sequences_no_repeats_hv1_ungapped.fa",
        hv2_ungapped="results/mash/opa_sequences_no_repeats_hv2_ungapped.fa"
    output:
        sv_sketch="results/mash/opa_sequences_no_repeats_sv_ungapped.fa.msh",
        hv1_sketch="results/mash/opa_sequences_no_repeats_hv1_ungapped.fa.msh",
        hv2_sketch="results/mash/opa_sequences_no_repeats_hv2_ungapped.fa.msh",
        sv_dist="results/mash/sv_mash_distance.tab",
        hv1_dist="results/mash/hv1_mash_distance.tab",
        hv2_dist="results/mash/hv2_mash_distance.tab"
    conda:
        "../envs/mash.yml"
    shell:
        """
            # Sketching
            mash sketch -i -k 6 {input.sv_ungapped}
            mash sketch -i -k 7 {input.hv1_ungapped}
            mash sketch -i -k 7 {input.hv2_ungapped}

            # Calculate distance matrix
            mash dist {output.sv_sketch} {output.sv_sketch} > {output.sv_dist}
            mash dist {output.hv1_sketch} {output.hv1_sketch} > {output.hv1_dist}
            mash dist {output.hv2_sketch} {output.hv2_sketch} > {output.hv2_dist}
        """

rule mcl:
    input:
        "results/mash/{region}_mash_distance.tab"
    output:
        abc="results/mcl/{region}_mash_distance.abc",
        mci="results/mcl/{region}.mci",
        tab="results/mcl/{region}.tab",
        inflation_mci=expand("results/mcl/out.{{region}}.mci.I{inflation}", inflation=INFLATIONS),
        dump=expand("results/mcl/dump.{{region}}.nci.I{inflation}", inflation=INFLATIONS),
        info="results/mcl/{region}.info",
        meet="results/mcl/{region}.meet",
        meet_dist="results/mcl/{region}.meet_dist"
    conda:
        "../envs/mcl.yml"
    shell:
        """
            mkdir -p results/mcl

            # First reformat the p-values from MASH
            cut -f 1,2,4 {input} > {output.abc}

            # Next write it into the MCL format (and transform the p-values to be -log10(p-value))
            mcxload -abc {output.abc} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {output.mci} -write-tab {output.tab}

            # Run MCL
            mcl {output.mci} -I 1.4 -o results/mcl/out.{wildcards.region}.mci.I14
            mcl {output.mci} -I 2 -o results/mcl/out.{wildcards.region}.mci.I20
            mcl {output.mci} -I 4 -o results/mcl/out.{wildcards.region}.mci.I40
            mcl {output.mci} -I 6 -o results/mcl/out.{wildcards.region}.mci.I60
            mcl {output.mci} -I 8 -o results/mcl/out.{wildcards.region}.mci.I80
            mcl {output.mci} -I 10 -o results/mcl/out.{wildcards.region}.mci.I100
            mcl {output.mci} -I 12 -o results/mcl/out.{wildcards.region}.mci.I120
            mcl {output.mci} -I 14 -o results/mcl/out.{wildcards.region}.mci.I140

            # Get the labels of the clusters
            mcxdump -icl results/mcl/out.{wildcards.region}.mci.I14 -tabr {output.tab} -o results/mcl/dump.{wildcards.region}.nci.I14
            mcxdump -icl results/mcl/out.{wildcards.region}.mci.I20 -tabr {output.tab} -o results/mcl/dump.{wildcards.region}.nci.I20
            mcxdump -icl results/mcl/out.{wildcards.region}.mci.I40 -tabr {output.tab} -o results/mcl/dump.{wildcards.region}.nci.I40
            mcxdump -icl results/mcl/out.{wildcards.region}.mci.I60 -tabr {output.tab} -o results/mcl/dump.{wildcards.region}.nci.I60
            mcxdump -icl results/mcl/out.{wildcards.region}.mci.I80 -tabr {output.tab} -o results/mcl/dump.{wildcards.region}.nci.I80
            mcxdump -icl results/mcl/out.{wildcards.region}.mci.I100 -tabr {output.tab} -o results/mcl/dump.{wildcards.region}.nci.I100
            mcxdump -icl results/mcl/out.{wildcards.region}.mci.I120 -tabr {output.tab} -o results/mcl/dump.{wildcards.region}.nci.I120
            mcxdump -icl results/mcl/out.{wildcards.region}.mci.I140 -tabr {output.tab} -o results/mcl/dump.{wildcards.region}.nci.I140

            # Calculate the clustering info
            clm info -o {output.info} {output.mci} {output.inflation_mci}

            # Calcualte the cluster distances (between different clusters created with different parameters as well as the "meet", the intersection of all the different clusters)
            clm meet -o {output.meet} {output.inflation_mci}
            clm dist -o {output.meet_dist} {output.meet} {output.inflation_mci}
        """

rule mcl_summary:
    input:
        info="results/mcl/{region}.info",
        meet_dist="results/mcl/{region}.meet_dist"
    output:
        summary="results/mcl/{region}_summary.csv"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            workflow/scripts/summarize_mcl_clusters.py {input.info} {input.meet_dist} {output.summary}
        """

rule mcl_plot:
    input:
        summary="results/mcl/{region}_summary.csv"
    output:
        multiext("figures/clusters/mcl_cluster_summary_{region}", ".png", ".pdf")
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            mkdir -p figures/clusters
            workflow/scripts/plot_mcl_cluster_summary.py {input.summary} {output} {wildcards.region}
        """

rule get_mcl_clusters:
    input:
        "results/mcl/dump.sv.nci.I80",
        "results/mcl/dump.hv1.nci.I80",
        "results/mcl/dump.hv2.nci.I140"
    output:
        "results/mcl/mcl_clusters.csv"
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
            workflow/scripts/postprocess_mcl_clusters.py {input} {output}
        """
