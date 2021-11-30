# Example code for how to use the USC HPC to run R/JAGS in parallel  

## Background  
This repository provides an example of running R and JAGS in parallel on USC's Discovery high performance computing cluster (https://www.carc.usc.edu/user-information/user-guides/hpc-basics/getting-started-discovery). The USC HPC uses the slurm workload manager to allocate resources and submit jobs. This code was informed by the materials available on the USCHPC github page (https://github.com/uschpc/workshop-r-hpc), which has simple examples of how to submit jobs to the USC HPC using slurm. 

For this specific analysis, we use MPI and pbdMPI (https://cran.r-project.org/web/packages/pbdMPI/index.html) to parallelize the R/rjags code. This code is designed to minimize the run time of embarrassingly parallel (or embarrassingly parallel) computations. Specifically, the slurm job script spreads the computation across 8 nodes and 16 tasks per node. This results in 128 connections with 128 CPU cores, which is the maximum number of connections permissible in R (without rebuilding R: https://stat.ethz.ch/pipermail/r-sig-hpc/2012-May/001373.html). As a note, the R package slurmR, in conjunction with the future package or in conjunction with the doParallel package, does not seem to work for this analysis. 

The resulting output is written to a csv file, where each CPU core writes results into a group of rows, and results of each iteration are written into new columns. To give an example of what the output looks like, let’s say you’re are running a function which results in a dataframe with 5 columns and 10 rows. If you run this function 12,800 times, and this computation is spread across 128 cores, this means that each core will perform 100 iterations of the function. The output of this function will then be a csv file with 128*x*10 rows, and 100*x*10 columns.


## Running analysis   

The basics of running this analysis are as follows: 

1) Set up the HPC cluster to run the analysis, using the "Shell commands for setting up HPC prior to analysis.sh". This includes installing R2jags.  
2) Upload slurm_job_script.job and r script to the HPC. The easiest way to do this is to use the CARC on-demand interface (https://www.carc.usc.edu/user-information/user-guides/hpc-basics/getting-started-ondemand). 
3) Use the code from "Shell commands for HPC analysis.sh" to submit job. 
4) Troubleshoot errors generated during stages 1-3. Rinse and repeat. 
