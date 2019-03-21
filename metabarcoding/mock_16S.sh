#!/bin/bash
#$ -q normal.q
#$ -N mock_16S
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

echo "Activation du module conda"
module purge
module load system/conda/5.1.0

echo "Activation de l'environnement conda"
conda activate qiime2-2018.11

MOCK_S=/homedir/galati/mock/analysis/16S/single/
MOCK_P=/homedir/galati/mock/analysis/16S/pair/

for i in $RAW_16S*_L001_R1_001.fastq.gz
do
    prefix_16S=$(basename $i _L001_R1_001.fastq.gz) # print which sample is being processed
ls ${MOCK_S}${prefix_16S}*_R1*.fastq* ${MOCK_S}${prefix_16S}*_R2*.fastq*

    echo "###Trim galore Sample" $prefix_16S "Start###"
echo "Import des fichiers fastq.gz"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${MOCK_S}${prefix_16S}* \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path 16S_demux-paired-end.qza
    echo "###Trim galore Sample" $prefix_16S "DONE###"
    
done
