# 
#$ -cwd 
#$ -j y 
#$ -pe orte 12
#$ -S /bin/bash 
#
# Name the job #$ -N RScript #
export PATH=/usr/local/Revo_rh6/bin:$PATH
Rscript ../R/main_sim_detail_q_cluster.R
echo "" 
echo "Done at " `date`