# 
#$ -cwd 
#$ -j y 
#$ -pe orte 24
#$ -S /bin/bash 
#
# Name the job #$ -N RScript #
export PATH=/usr/local/Revo_rh6/bin:$PATH
Rscript ../R/MainExpRealCompareLR.R
echo ""
echo "Done at " `date`