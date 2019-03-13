#!/bin/bash
#$ -q normal.q
#$ -N trimgalore
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

echo "Activation du module conda"
module load system/conda/5.1.0

echo "Activation de l'environnement conda"
conda activate qiime2-2018.1

echo "Suppression des adapters ITS"
RAW_ITS=/homedir/galati/data/ITS/
TRIM_ITS=/homedir/galati/data/ITS_adapter/

for i in $RAW_ITS*_L001_R1_001.fastq.gz
do
    prefix_ITS=$(basename $i _L001_R1_001.fastq.gz) # print which sample is being processed
ls ${RAW_ITS}${prefix_ITS}*_R1*.fastq* ${RAW_ITS}${prefix_ITS}*_R2*.fastq*

    echo "###Trim galore Sample" $prefix_ITS "Start###"
trim_galore --paired -q 0 --nextera --length 0 $RAW_ITS${prefix_ITS}*_R1*.fastq* $RAW_ITS${prefix_ITS}*_R2*.fastq* -o $TRIM_ITS
    echo "###Trim galore Sample" $prefix_ITS "DONE###"

done

echo "Suppression des adapters 16S"
RAW_16S=/homedir/galati/data/16S/
TRIM_16S=/homedir/galati/data/16S_adapter/

for i in $RAW_16S*_L001_R1_001.fastq.gz
do
    prefix_16S=$(basename $i _L001_R1_001.fastq.gz) # print which sample is being processed
ls ${RAW_16S}${prefix_16S}*_R1*.fastq* ${RAW_16S}${prefix_16S}*_R2*.fastq*

    echo "###Trim galore Sample" $prefix_16S "Start###"
trim_galore --paired -q 0 --nextera --length 0 $RAW_16S${prefix_16S}*_R1*.fastq* $RAW_16S${prefix_16S}*_R2*.fastq* -o $TRIM_16S
    echo "###Trim galore Sample" $prefix_16S "DONE###"
    
done

for f in ${TRIM_ITS}*.fq.gz ; do mv $f $(echo $f | cut -d "_" -f 1,2,3,4,5).fastq.gz ; done ;
for f in ${TRIM_16S}*.fq.gz ; do mv $f $(echo $f | cut -d "_" -f 1,2,3,4,5).fastq.gz ; done ;

echo "Import des fichiers fastq.gz"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /homedir/galati/data/16S_adapter \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /homedir/galati/data/16S_demux-paired-end.qza

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /homedir/galati/data/ITS_adapter \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path /homedir/galati/data/ITS_demux-paired-end.qza


echo "Suppression des primers 16S"
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux-paired-end.qza \
        --p-cores 1 \
        --p-front-f NNNNCCTACGGGNGGCWGCAG \
        --p-front-r NNNNGACTACHVGGGTATCTAATCC \
        --p-error-rate 0.1 \
        --o-trimmed-sequences demux-paired-end-trimmed.qza \
        
echo "Suppression des primers ITS"
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux-paired-end.qza \
        --p-cores 1 \
        --p-front-f NNNNGTGAATCATCGAATCTTTGAA \
        --p-front-r NNNNTCCTCCGCTTATTGATATGC \
        --p-error-rate 0.1 \
        --o-trimmed-sequences demux-paired-end-trimmed.qza \
        
