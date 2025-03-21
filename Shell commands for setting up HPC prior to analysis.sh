
#############################################################
# Commands for troubleshooting and testing ------------------

# Log in
ssh jagoodri@endeavour.usc.edu



# TO SOLVE ERROR WHERE JAGS IS NOT RECOGNIZED ON WORKER NODES:
# (this only needs to be run once)
module load gcc/13.3.0 patchelf
patchelf --set-rpath /spack/apps2/linux-centos7-x86_64/gcc-13.3.0/jags-4.3.0-d373ytkqeamusbww7n2qjxyfwgikw2w5/lib /home1/jagoodri/R/x86_64-pc-linux-gnu-library/4.1/rjags/libs/rjags.so



# To install R2jags
module spider jags
module load usc r
module load jags
module load gcc/13.3.0
module load openblas/0.3.28


export PKG_CONFIG_PATH=/packages/jags/4.3.0/lib/pkgconfig
pkg-config ––modversion jags
R
install.packages("rjags", configure.args="––enable-rpath")
install.packages("R2jags", configure.args="––enable-rpath")
