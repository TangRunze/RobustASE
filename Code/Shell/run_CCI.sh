# 
#$ -cwd 
#$ -j y 
#$ -pe orte 24
#$ -S /bin/bash 
#
# Name the job #$ -N RScript #
export PATH=/usr/local/Revo_rh6/bin:$PATH
Rscript ../R/MainExpCCI2Version.R 
echo "" 
echo "Done at " `date`