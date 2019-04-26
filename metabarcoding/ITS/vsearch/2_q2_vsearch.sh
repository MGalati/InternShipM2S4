#!/bin/bash

#$ -q bigmem.q
#$ -N vsearch_ITS
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

#https://forum.qiime2.org/t/analyzing-paired-end-reads-in-qiime-2/2021

module purge
module load system/conda/5.1.0
source activate qiime2-2018.11

# JOB BEGIN

IN=/homedir/galati/data/metab/ITS

RUN1=PRIM
RUN2=ITS_mock24
"""
echo 'Import'
for seqs in ${RUN1} ${RUN2}
do
qiime tools import --type SampleData[PairedEndSequencesWithQuality] \
                   --input-path ${IN}/${seqs} \
                   --output-path ${IN}/${seqs}_reads.qza \
                   --input-format CasavaOneEightSingleLanePerSampleDirFmt 

#Check this artifact to make sure that QIIME now recognizes your data
qiime tools peek ${IN}/${seqs}_reads.qza 

echo 'Initial sequence quality control'
qiime demux summarize \
  --i-data ${IN}/${seqs}_reads.qza  \
  --o-visualization ${IN}/${seqs}_reads.qzv  \
  --verbose
done
"""
echo 'Début de boucle vsearch'
for seqs in ${RUN1} ${RUN2}
do

mkdir vsearch_output
mkdir vsearch_output_${seqs}

echo 'Quality filter'
qiime quality-filter q-score \
 --i-demux ${IN}/${seqs}_reads.qza \
 --o-filtered-sequences vsearch_output_${seqs}/${seqs}_demux-filtered.qza \
 --o-filter-stats vsearch_output_${seqs}/${seqs}_demux-filter-stats.qza

##Parameter 'demultiplexed_seqs' received an argument of type SampleData[SequencesWithQuality]. An argument of subtype SampleData[PairedEndSequencesWithQuality] is required.
##echo 'Joining reads'
##qiime vsearch join-pairs \
##  --i-demultiplexed-seqs vsearch_output_${seqs}/${seqs}_demux-filtered.qza \
##  --o-joined-sequences ${IN}/${seqs}_demux-joined.qza

echo 'Summurize joined data with read quality'
qiime demux summarize \
  --i-data ${IN}/${seqs}_demux-joined.qza \
  --o-visualization ${IN}/${seqs}_demux-joined.qzv
  
echo 'Dereplicate sequences'
qiime vsearch dereplicate-sequences \
  --i-sequences ${IN}/${seqs}_demux-joined.qza \
  --o-dereplicated-table ${IN}/${seqs}_table.qza \
  --o-dereplicated-sequences ${IN}/${seqs}_rep-seqs.qza

echo 'De novo clustering'
qiime vsearch cluster-features-de-novo \
  --i-table ${IN}/${seqs}_table.qza \
  --i-sequences ${IN}/${seqs}_rep-seqs.qza \
  --o-clustered-table vsearch_output_${seqs}/${seqs}_table-dn-99.qza \
  --o-clustered-sequences vsearch_output_${seqs}/${seqs}_rep-seqs-dn-99.qza \
  --p-perc-identity 0.99

echo 'summarize your filtered/ASV table data'
qiime tools export --input-path vsearch_output_${seqs}/${seqs}_demux-filter-stats.qza --output-path vsearch_output_${seqs}/${seqs}

qiime feature-table summarize --i-table ${IN}/${seqs}_table.qza --o-visualization ${IN}/${seqs}_table_summary.qzv --verbose

done
echo 'Fin de boucle'
echo 'Merging denoised data'

echo 'ASV table'
qiime feature-table merge \
  --i-tables vsearch_output_${RUN1}/${RUN1}_table-dn-99.qza \
  --i-tables vsearch_output_${RUN2}/${RUN2}_table-dn-99.qza \
  --o-merged-table vsearch_output/table.qza

echo 'summarize'
qiime feature-table summarize \
  --i-table vsearch_output/table.qza \
  --o-visualization vsearch_output/table.qzv 

echo 'Representative sequences'
qiime feature-table merge-seqs \
  --i-data vsearch_output_${RUN1}/${RUN1}_rep-seqs-dn-99.qza \
  --i-data vsearch_output_${RUN2}/${RUN2}_rep-seqs-dn-99.qza \
  --o-merged-data vsearch_output/representative_sequences.qza

qiime feature-table tabulate-seqs \
  --i-data dada2_output/representative_sequences.qza\
  --o-visualization dada2_output/representative_sequences.qzv

mkdir phylogeny

echo 'Phylogenie'
qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads ${NSLOTS} \
  --i-sequences vsearch_output/representative_sequences.qza \
  --o-alignment phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree phylogeny/unrooted-tree.qza \
  --o-rooted-tree phylogeny/rooted-tree.qza \
  --verbose

### Assign Taxonomy
# loop to test various taxonomic database - pour toi laisser juste silva123 - 
# https://www.dropbox.com/s/5tckx2vhrmf3flp/silva-132-99-nb-classifier.qza?dl=0

mkdir taxonomy
echo 'Classifier'
qiime feature-classifier classify-sklearn \
  --i-classifier /homedir/galati/data/metab/ITS/classifier/unite-ver7-dynamic-classifier-01.12.2017.qza \
  --i-reads vsearch_output/representative_sequences.qza \
  --o-classification taxonomy/ITS_taxonomy.qza \
  --p-n-jobs ${NSLOTS} \
  --verbose

qiime metadata tabulate \
  --m-input-file taxonomy/ITS_taxonomy.qza \
  --o-visualization taxonomy/ITS_taxonomy.qzv

echo 'necessite metadata'
 qiime taxa barplot \
  --i-table vsearch_output/table.qza \
  --i-taxonomy taxonomy/ITS_taxonomy.qza \
  --o-visualization taxonomy/ITS_taxa-bar-plots.qzv \
  --m-metadata-file /homedir/galati/data/metab/16S/metadata/metadata.tsv 

qiime tools export --input-path taxonomy/ITS_taxonomy.qza --output-path taxonomy
mv taxonomy/taxonomy.tsv taxonomy/ITS_taxonomy.tsv


### Exporting and modifying BIOM tables

#Creating a TSV BIOM table
qiime tools export --input-path vsearch_output/table.qza --output-path export
biom convert -i export/feature-table.biom -o export/ASV-table.biom.tsv --to-tsv

cp export/ASV-table.biom.tsv export/feature-table.biom.tsv

#https://askubuntu.com/questions/20414/find-and-replace-text-within-a-file-using-commands
# ASV table for phyloseq

sed -i "s/#OTU ID/OTUID/g" export/ASV-table.biom.tsv
sed -i "1d" export/ASV-table.biom.tsv

sed -i "s/#OTU ID/#OTUID/g" export/feature-table.biom.tsv

#Export Taxonomy
qiime tools export --input-path /homedir/galati/data/metab/ITS/classifier/unite-ver7-dynamic-classifier-01.12.2017.qza --output-path export

biom add-metadata -i export/ASV-table.biom.tsv  -o export/ASV-table-unite-ver7-taxonomy.biom \
  --observation-metadata-fp export/unite-ver7_taxonomy.tsv \
  --sc-separated taxonomy
biom convert -i export/ASV-table-unite-ver7-taxonomy.biom -o export/ASV-table-unite-ver7-taxonomy.biom.tsv --to-tsv

#Export ASV seqs
qiime tools export --input-path vsearch_output/representative_sequences.qza --output-path export

#Export Tree
qiime tools export --input-path phylogeny/unrooted-tree.qza --output-path export
mv export/tree.nwk export/unrooted-tree.nwk
qiime tools export --input-path phylogeny/rooted-tree.qza --output-path export
mv export/tree.nwk export/rooted-tree.nwk


# For phyloseq you will need :
#ls export/feature-table.biom.tsv export/taxonomy.tsv export/unrooted-tree.nwk export/rooted-tree.nwk

zip export/export.zip export/* vsearch_outpu*/*qzv taxonomy/*.qzv


# JOB END
date

exit 0
