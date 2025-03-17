rule get_inputs_for_beast:
    input:
        metadata='input_data/lab_strains_hybrid_genomes_metadata.csv',
        alignment='results/gubbins/complete_genome_pseudogenomes.filtered_polymorphic_sites.fasta'
    output:
        alignment_with_dates='results/beast/complete_genomes_recomb_masked_with_dates.fasta',
        dates='results/beast/dates.dat'
    conda:
        '../envs/biopython_plotting.yml'
    shell:
        """
            mkdir -p results/beast
            workflow/scripts/get_inputs_for_beast.py {input.metadata} {input.alignment} {output.alignment_with_dates} {output.dates}
        """