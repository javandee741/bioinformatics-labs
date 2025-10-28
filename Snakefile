# Snakefile (c маячком)
rule all:
    input:
        "results/amino_stats.tsv",
        "results/plots.done"

rule amino_stats:
    input:
        fasta="files/protein.fasta"
    output:
        tsv="results/amino_stats.tsv",
        plots_done="results/plots.done"
    params:
        plots_dir="results"
    shell:
        """
        python scripts/amino_counter.py \
            -i {input.fasta} \
            -o {output.tsv} \
            --plots-dir {params.plots_dir}
        touch {output.plots_done}
        """
