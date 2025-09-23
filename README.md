This is a piece of modular code for defining reaction networks, placing enzymes within a spherical object and computing the steady state spatial distribution of the reactants, intermediates, and products.

NOT SURE ABOUT THIS. PRELIMINARY!

Best to use conda:
1. make a new environment named snakemake-runner
    ```bash
    conda create -n snakemake-runner -c bioconda snakemake
    ```
2. activate the new environment
    ```bash
    conda activate snakemake-runner
    ```
3. run the snakemake files
    ```bash
    run python -m snakemake -s Snakefile.smk check_case_input_validity --config case_number=0 --cores 1 --use-conda
    ```