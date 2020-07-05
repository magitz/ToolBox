#!/bin/bash
#SBATCH --qos=soltis-b                 # Your qos.
#SBATCH --account=soltis               # Your account or group.
#SBATCH --mail-user=USER@ufl.edu       # Change to your email or omit
#SBATCH --time=48:00:00                # Change to max time estimated for longest job
###############################################################################
####### Don't change the info in the %% on the lines below, 
####### raxml_job_generator.py will change them for you.
###############################################################################
#SBATCH --cpus-per-task=%%THREADS%%    # Will be replaced by raxml_job_generator.py
#SBATCH --nodes=1                      # raxml_job_generator.py does not do MPI jobs, leave at 1
#SBATCH --ntasks=1                     # raxml_job_generator.py does not do MPI jobs, leave at 1
#SBATCH --mem=%%MEMORY%%MB             # Will be replaced by raxml_job_generator.py
#SBATCH --job-name=%%JOBNAME%%         # Will be replaced by raxml_job_generator.py
#SBATCH --output=%%JOB_OUT%%           # Will be replaced by raxml_job_generator.py
#SBATCH --error=%%JOB_ERROR%%          # Will be replaced by raxml_job_generator.py

# %%INFO%%

#Run raxml
module load raxml-ng

raxml-ng --search --msa %%ALIGNMENT%% --threads %%THREADS%% --model %%MODEL%% --prefix %%PREFIX%%
