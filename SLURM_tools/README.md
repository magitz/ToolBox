# SLURM Tools

I have just started adding scripts to this directory, but expect to add more as time goes on.

## raxml_job_generator.py

While job arrays can be handy, they suffer from the fact that all tasks in the array have the same resource request. With multi-gene datasets, this isn't always the most efficient method as you need to request the resources for the largest job and thus smaller jobs will waste resources (too much memory, too many cores).

[RAxML-ng](https://github.com/amkozlov/raxml-ng) has a handy `--parse` option which estimates the RAM and CPU threads needed to analyze a dataset. The `raxml_job_generator.py` script automates the task of running `raxml-ng --parse` on all of the `.phy` files in a directory, parsing the RAM and CPU estimates, using a template SLURM script to substitute customized values, and optionally, submitting the job scripts. Thus, the convenience of job arrays can be achieved, but with the efficiency of per-dataset resource requests.

### Usage

```bash
raxml_job_generator.py [-h] -d DIRECTORY -o OUTPUT_DIR [-l LOG_DIR]
                              [-j JOB_NAME_PREFIX] [-m MODEL] [-t TEMPLATE]
                              [-c MAX_CPU] [-b MEM_BUFFER] [-r RAXML] [-n]
                              [-s] [-v] [--version]
```

Argument   | Description
-------|------------
-h/--help | Print help message
**-d/--directory** | Directory where alignment files are located. These should end in `.phy`. *This argument is required.*
**-o/--output_dir** | Relative path to output directory. This is where the submit scripts will be saved and submitted from, if `--submit` option is used. *This argument is required*.
-l/--log_dir | Name for the directory to put SLURM job output and error files. Will be nested within the `--output_dir`. Default: `logs`.
-j/--job_name_prefix | Name for RAxML-ng `--prefix` option, this plus the alignment filename will be used. Default: use only alignment name.
-m/--model | Evolutionary model to use for the RAxML-ng `--model` option. Default: `GTR+G`.
-t/--template | Path to SLURM template submit script. This script will be used to generate all of the customized scripts using the text marked with %%OPTION%%. Default: [RAxML_Gene_Tree.sbatch](RAxML_Gene_Tree.sbatch) script in this repository.
-c/--max_cpu | The maximum number of CPUs (cores) to allocate per job. If the RAxML-ng estimate is higher than this, the job script will limit to this value. Default: `32`.
-b/--mem_buffer | Memory buffer to add to job memory requests over what RAxML-ng suggests. This is a percentage, e.g. 0.15 is 15%. **Note** a minimum of 100MB will always be requested. Default: `0.15`.
-r/--raxml | path to raxml-ng executable. Default: `raxml-ng` (either same directory or in $PATH)
-s/--submit | Submit the jobs as they are generated. If not passed, job scripts are saved, but not submitted.
-v/--verbose | Verbose output.

This script was developed for running HiPerGator at the University of Florida, but should generally work for most SLURM implementations. You may need to modify the template script more heavily for other cluster setups. If you modify the script, take care not to change the customization text indicated with %%OPTION%% marks.

Thanks to [Heather Kates](https://github.com/HeatherKates) for the suggestion for this script and her original bash implementation.
