# Log in to HPC 
# Performed analysis on the epyc-64	partition on Discovery

# Log ig
ssh jagoodri@discovery.usc.edu

# Change Working Directory
cd /project/dconti_624/Users/jagoodri/sol


# Run analysis
dos2unix example_slurm_job_script.job
sbatch example_slurm_job_script.job

exit
