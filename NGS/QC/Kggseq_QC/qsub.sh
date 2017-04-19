#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l mem=16gb
#PBS -l walltime=12:00:00
#PBS -m a
#PBS -q medium
#PBS -N KGGseq_QC


cd $PBS_O_WORKDIR

Path=$(pwd)
txt=$Path/Kggseq_basicQC.txt
kggseq=/psychipc01/disk2/software/KGGseq/20170401/kggseq10/kggseq_20170401.jar

/psychipc01/disk2/software/java-8/jdk1.8.0_111/bin/java -jar $kggseq $txt
