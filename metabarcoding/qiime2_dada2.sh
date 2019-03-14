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

echo 'Renommage des fichiers au format Casava et à l'extension .fastq.gz"
for f in ${TRIM_ITS}*.fq.gz ; do mv $f $(echo $f | cut -d "_" -f 1,2,3,4,5).fastq.gz ; done ;
for f in ${TRIM_16S}*.fq.gz ; do mv $f $(echo $f | cut -d "_" -f 1,2,3,4,5).fastq.gz ; done ;

echo "Suppresison des fichiers de rapport"
rm ${TRIM_ITS}*.txt
rm ${TRIM_ITS}*.txt

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
        --i-demultiplexed-sequences 16S_demux-paired-end.qza \
        --p-cores 1 \
        --p-front-f NNNNCCTACGGGNGGCWGCAG \
        --p-front-r NNNNGACTACHVGGGTATCTAATCC \
        --p-error-rate 0.1 \
        --o-trimmed-sequences 16S_demux-paired-end-trimmed.qza \
        
echo "Suppression des primers ITS"
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences ITS_demux-paired-end.qza \
        --p-cores 1 \
        --p-front-f NNNNGTGAATCATCGAATCTTTGAA \
        --p-front-r NNNNTCCTCCGCTTATTGATATGC \
        --p-error-rate 0.1 \
        --o-trimmed-sequences ITS_demux-paired-end-trimmed.qza \

echo "Fichiers de visualisation des données"
qiime demux summarize --i-data ITS_demux-paired-end.qza --o-visualization ITS_demux.qzv
qiime demux summarize --i-data 16S_demux-paired-end.qza --o-visualization 16S_demux.qzv
scp galati@cc2-login.cirad.fr:/homedir/galati/data/ITS_demux.qzv /home/galati/Téléchargements/
scp galati@cc2-login.cirad.fr:/homedir/galati/data/16S_demux.qzv /home/galati/Téléchargements/

echo "Denoising 16S"
qiime dada2 denoise-paired --i-demultiplexed-seqs 16S_demux-paired-end.qza \
                           --p-trunc-len-f 0 \
                           --p-trunc-len-r 207 \ #Quality Score decreases from 240pb for the reverse reads
                           --p-max-ee 2.0 \ #default value : all the reads with number of exepcted errors higher than 2.0 will be discarded
                           --p-trunc-q 10 \ #reads are truncated at the first instance of a quality score less than or equal to 10
                           --p-n-reads-learn 1000000 \ #default value : it's the number of read to use during the training of error model 
                           --p-n-threads ${NSLOTS} \ #default value. When =0, it uses all the available cores
                           --p-chimera-method consensus\ #default value : chimeras are detected in samples individually, and sequences found chimeric in a sufficient fraction of samples are removed
                           --o-representative-sequences 16S_rep-seq-dada2.qza \
                           --o-table 16S_table-dada2.qza \
                           --o-denoising-stats 16S_stats-dada2.qza \
                           --verbose \

echo "Denoising ITS"
qiime dada2 denoise-paired --i-demultiplexed-seqs ITS_demux-paired-end.qza \
                           --p-trunc-len-f 0 \
                           --p-trunc-len-r 220 \ #Quality Score decreases from 240pb for the reverse reads
                           --p-max-ee 2.0 \ #default value : all the reads with number of exepcted errors higher than 2.0 will be discarded
                           --p-trunc-q 10 \ #reads are truncated at the first instance of a quality score less than or equal to 10
                           --p-n-reads-learn 1000000 \ #default value : it's the number of read to use during the training of error model 
                           --p-n-threads ${NSLOTS} \ #default value. When =0, it uses all the available cores
                           --p-chimera-method consensus\ #default value : chimeras are detected in samples individually, and sequences found chimeric in a sufficient fraction of samples are removed
                           --o-representative-sequences ITS_rep-seq-dada2.qza \
                           --o-table ITS_table-dada2.qza \
                           --o-denoising-stats ITS_stats-dada2.qza \
                           --verbose \
        
echo "Phylogénie 16S"
qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads 0 \
  --i-sequences 16S_rep-seqs-dada2.qza \
  --o-alignment phylogeny/16S_aligned-rep-seqs.qza \
  --o-masked-alignment phylogeny/16S_masked-aligned-rep-seqs.qza \
  --o-tree phylogeny/16S_unrooted-tree.qza \
  --o-rooted-tree phylogeny/16S_rooted-tree.qza \
  --verbose \
  
echo "Phylogénie ITS"  
qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads 0 \
  --i-sequences ITS_rep-seqs-dada2.qza \
  --o-alignment phylogeny/ITS_aligned-rep-seqs.qza \
  --o-masked-alignment phylogeny/ITS_masked-aligned-rep-seqs.qza \
  --o-tree phylogeny/ITS_unrooted-tree.qza \
  --o-rooted-tree phylogeny/ITS_rooted-tree.qza \
  --verbose \
