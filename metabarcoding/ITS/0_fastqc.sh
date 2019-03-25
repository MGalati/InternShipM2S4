#!/bin/bash

#$ -q normal.q
#$ -N fastqc_ITS
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 6
#$ -l mem_free=6G
#$ -V
#$ -cwd
#$ -V


module purge
module load system/conda/5.1.0
# conda create --name fastqc_multiqc
conda activate fastqc_multiqc

echo "Début des analyses FastQC/MultiQC"

# Les fichiers doivent respecter le format Qiime2 :
# par exemple : L2S357_15_L001_R1_001.fastq.gz

MOCK=/homedir/galati/mock/analysis/ITS/pair/ITS_mock24
RAW_ITS=/homedir/galati/data/ITS/

for f in ${MOCK} ${RAW_ITS}
do
mkdir ${f}/qc/
fastqc -t $NSLOTS ${f}/*.fastq.gz -o ${f}/qc/
multiqc ${f}/qc/ -o ${f}/qc/
done

echo 'Fin des analyses'

date
 
exit 0
