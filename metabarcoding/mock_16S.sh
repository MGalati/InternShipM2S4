#!/bin/bash
#$ -q normal.q
#$ -N mock_16S
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

MOCK_S=/homedir/galati/mock/analysis/16S/single/
MOCK_P=/homedir/galati/mock/analysis/16S/pair/

for i in ${RAW_16S}*_L001_R1_001.fastq.gz
do
    prefix_16S=$(basename $i _L001_R1_001.fastq.gz) # print which sample is being processed
ls ${MOCK_S}${prefix_16S}*_R1*.fastq*

echo "Import des fichiers fastq.gz"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ${MOCK_S}${i}* --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path 16S_${i}_demux-paired-end.qza
    
done

# JOB END
date

exit 0
