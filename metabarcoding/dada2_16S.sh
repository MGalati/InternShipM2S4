#!/bin/bash
#$ -q normal.q
#$ -N dada2_16S
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
for f in ${TRIM_16S}*.fq.gz ; do mv $f $(echo $f | cut -d "_" -f 1,2,3,4,5).fastq.gz ; done ;

echo "Suppresison des fichiers de rapport"
rm ${TRIM_16S}*.txt

echo "Import des fichiers fastq.gz"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /homedir/galati/data/16S_adapter \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /homedir/galati/data/16S_demux-paired-end.qza
  
  echo "Suppression des primers 16S"
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences 16S_demux-paired-end.qza \
        --p-cores 1 \
        --p-front-f NNNNCCTACGGGNGGCWGCAG \
        --p-front-r NNNNGACTACHVGGGTATCTAATCC \
        --p-error-rate 0.1 \
        --o-trimmed-sequences 16S_demux-paired-end-trimmed.qza \

echo "Denoising 16S"
qiime dada2 denoise-paired --i-demultiplexed-seqs 16S_demux-paired-end.qza \
                           --p-trunc-len-f 0 \
                           --p-trunc-len-r 207 \ #Quality Score decreases from 207pb for the reverse reads
                           --p-max-ee 2.0 \ #default value : all the reads with number of exepcted errors higher than 2.0 will be discarded
                           --p-trunc-q 10 \ #reads are truncated at the first instance of a quality score less than or equal to 10
                           --p-n-reads-learn 1000000 \ #default value : it's the number of read to use during the training of error model 
                           --p-n-threads ${NSLOTS} \ #default value. When =0, it uses all the available cores
                           --p-chimera-method consensus\ #default value : chimeras are detected in samples individually, and sequences found chimeric in a sufficient fraction of samples are removed
                           --o-representative-sequences 16S_rep-seq-dada2.qza \
                           --o-table 16S_table-dada2.qza \
                           --o-denoising-stats 16S_stats-dada2.qza \
                           --verbose \
                           
mkdir /homedir/galati/data/phylo         
DATA=/homedir/galati/data
PHYLO=/homedir/galati/data/phylo       

echo "Phylogénie 16S"
qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads 0 \
  --i-sequences ${DATA}/16S_rep-seq-dada2.qza \
  --o-alignment ${PHYLO}/16S_aligned-rep-seqs.qza \
  --o-masked-alignment ${PHYLO}/16S_masked-aligned-rep-seqs.qza \
  --o-tree ${PHYLO}/16S_unrooted-tree.qza \
  --o-rooted-tree ${PHYLO}/16S_rooted-tree.qza \
  --verbose \
  
mkdir /homedir/galati/data/taxonomy/  

echo "Taxonomie des séquences 16S"
qiime feature-classifier classify-sklearn \
  --i-classifier /homedir/galati/data/classifier/silva-132-99-nb-classifier.qza \
  --i-reads /homedir/galati/data/16S_rep-seq-dada2.qza  \
  --o-classification /homedir/galati/data/taxonomy/16S_taxonomy.qza \
  --p-n-jobs 15 \
  --verbose
  
echo "Tabulate 16S"
qiime metadata tabulate \
  --m-input-file /homedir/galati/data/taxonomy/16S_taxonomy.qza  \
  --o-visualization /homedir/galati/data/taxonomy/16S_taxonomy.qzv  

# JOB END
date

exit 0
