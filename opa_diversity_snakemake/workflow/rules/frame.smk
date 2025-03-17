rule plot_opa_frame:
    input:
        "results/opa_metadata_locus.csv"
    output:
        expand("figures/frame/frame_on_by_strain.{figure_format}", figure_format=FIGURE_FORMATS)
    conda:
        "../envs/python_minimal.yml"
    shell:
        """
        mkdir -p figures/frame
        workflow/scripts/plot_opa_frame.py {input} {output}
        """