#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l mem=16gb
#PBS -l walltime=12:00:00
#PBS -m a
#PBS -q medium
#PBS -N KGGseq_QC

Path=/home/jacobhsu/scripts/shared-scripts
txt=$Path/NGS/QC/Kggseq_QC/Kggseq_basicQC.txt
kggseq=/home/groups/pcsham/shared/kggseq_V1/20160816/kggseq20160924.jar

cd $Path/NGS/QC/Kggseq_QC &&
java -jar $kggseq $txt
