This is a piece of modular code for defining reaction networks, placing enzymes within a spherical object and computing the steady state spatial distribution of the reactants, intermediates, and products.

NOT SURE ABOUT THIS. PRELIMINARY!

Best to use conda:
1. make a new environment named snakemake-runner
    ```bash
    conda create -n snakemake-env -c conda-forge -c bioconda snakemake
    ```
2. activate the new environment
    ```bash
    conda activate snakemake-env
    ```
3. run the snakemake files
    ```bash
    run python -m snakemake -s Snakefile.smk check_case_input_validity --config case_number=0 --cores 1 --use-conda
    ```

In case the command prompt does not recognize conda (on Windows), follow the next steps:

1. Open Anaconda Prompt
2. Check Conda installed location: `where conda`
3. Open Advanced System Settings
4. Click on Environment Variables
5. Edit Path
6. Copy paste what was found in `where conda`

if it does not work, edit path with `C:\Users\YourUsername\Anaconda3\Scripts` or analogous

Important:

(Taken from https://github.com/snakemake/snakemake/issues/1619): Note that conda environments are only used with shell, script, notebook and the wrapper directive, not the run directive. The reason is that the run directive has access to the rest of the Snakefile (e.g. globally defined variables) and therefore must be executed in the same process as Snakemake itself.