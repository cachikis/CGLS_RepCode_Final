#!/bin/bash
#PBS -N RunMainCGLS
#PBS -l nodes=1:ppn=15,mem=128gb
#PBS -j oe
#PBS -o /home/cachikis/CGLStest2/CGLS_RepCode_Final/Logs/output.txt


# remember to change the output path so that is comports with
# your local file system

# Change "cd" statement to your own directory
cd "/home/cachikis/CGLStest2/CGLS_RepCode_Final/"
workingdirectory=/home/cachikis/CGLStest2/CGLS_RepCode_Final/


# execute program
Rscript Code/Data/data_filter_v4.R $workingdirectory
Rscript Code/Data/gomme1.R $workingdirectory
Rscript Code/Data/gomme_process_IA.R $workingdirectory

# Matlab can clobber it's temporary files if multiple instances are run at once
# Create a job specific temp directory to avoid this
mkdir -p ~/matlabtmp/$PBS_JOBID  
export MATLABWORKDIR=~/matlabtmp/$PBS_JOBID

cd Code/Main/ 
matlab -nodesktop -nosplash -r "run('main_cgls.m'); exit;"
cd ../../

# clean up
rm -rf ~/matlabtmp
