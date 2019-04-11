#!/bin/bash

#$ -q bigmem.q
#$ -N QC_16S
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 6
#$ -l mem_free=6G
#$ -V
#$ -cwd


module purge
module load system/conda/5.1.0
source activate fastqc_multiqc

echo "Début des analyses FastQC/MultiQC"

# Les fichiers doivent respecter le format Qiime2 :
# par exemple : L2S357_15_L001_R1_001.fastq.gz

MOCK=/homedir/galati/data/metab/16S/MOCK/Mock_S280
RAW_16S=/homedir/galati/data/metab/16S/RAW

for f in ${MOCK} ${RAW_16S}
do
mkdir ${f}/qc/
fastqc -t $NSLOTS ${f}/*.fastq.gz -o ${f}/qc/
multiqc ${f}/qc/ -o ${f}/qc/
done

echo 'Fin des analyses'

mkdir /homedir/galati/data/metab/16S/RAW_qc/
mv /homedir/galati/data/metab/16S/RAW/qc/* /homedir/galati/data/metab/16S/RAW_qc/

date
 
exit 0
