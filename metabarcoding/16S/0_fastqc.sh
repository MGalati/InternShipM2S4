#!/bin/bash

#$ -q normal.q
#$ -N fastqc_16S
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 6
#$ -l mem_free=6G
#$ -V
#$ -cwd
#$ -V


module purge
module load system/conda/5.1.0
# conda create --name fastqc_multiqc
source activate fastqc_multiqc

echo "Début des analyses FastQC/MultiQC"

# Les fichiers doivent respecter le format Qiime2 :
# par exemple : L2S357_15_L001_R1_001.fastq.gz

MOCK=/homedir/galati/mock/analysis/16S/pair/Mock_S280
RAW_16S=/homedir/galati/data/16S/

for f in ${MOCK} ${RAW_16S}
do
mkdir ${f}/qc/
fastqc -t $NSLOTS ${f}/*.fastq.gz -o ${f}/qc/
multiqc ${f}/qc/ -o ${f}/qc/
done

echo 'Fin des analyses'

date
 
exit 0
