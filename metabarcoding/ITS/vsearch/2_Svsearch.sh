#!/bin/bash

#$ -q short.q
#$ -N Svsearch_ITS
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 1
#$ -l mem_free=6G
#$ -V
#$ -cwd

#https://forum.qiime2.org/t/analyzing-paired-end-reads-in-qiime-2/2021

module purge
module load system/conda/5.1.0
source activate qiime2-2018.11

# JOB BEGIN

IN=/homedir/galati/data/metab/ITS

RUN1=SUB
RUN2=ITS_mock24

for seqs in ${RUN1} ${RUN2}
do
qiime tools import --type SampleData[PairedEndSequencesWithQuality] \
                   --input-path ${IN}/${seqs} \
                   --output-path ${IN}/${seqs}_reads.qza \
                   --input-format CasavaOneEightSingleLanePerSampleDirFmt 

#Check this artifact to make sure that QIIME now recognizes your data
qiime tools peek ${IN}/${seqs}_reads.qza 

### 'Initial' sequence quality control
qiime demux summarize \
  --i-data ${IN}/${seqs}_reads.qza  \
  --o-visualization ${IN}/${seqs}_reads.qzv  \
  --verbose
done

for seqs in ${RUN1} ${RUN2}
do

mkdir vsearch_output
mkdir vsearch_output_${seqs}

### Joining reads
qiime vsearch join-pairs \
  --i-demultiplexed-seqs ${IN}/${seqs}_reads.qza \
  --o-joined-sequences ${IN}/${seqs}_demux-joined.qza
  
### Summurize joined data with read quality
qiime demux summarize \
  --i-data ${IN}/${seqs}_demux-joined.qza \
  --o-visualization ${IN}/${seqs}_demux-joined.qzv

### Quality filter
qiime quality-filter q-score \
 --i-demux ${IN}/${seqs}_demux.qza \
 --o-filtered-sequences vsearch_output_${seqs}/${seqs}_demux-filtered.qza \
 --o-filter-stats vsearch_output_${seqs}/${seqs}_demux-filter-stats.qza

### Dereplicate sequences
qiime vsearch dereplicate-sequences \
  --i-sequences ${IN}/${seqs}_reads.qza \
  --o-dereplicated-table ${IN}/${seqs}_table.qza \
  --o-dereplicated-sequences ${IN}/${seqs}_rep-seqs.qza

### De novo clustering
qiime vsearch cluster-features-de-novo \
  --i-table ${IN}/${seqs}_table.qza \
  --i-sequences ${IN}/${seqs}_rep-seqs.qza \
  --o-clustered-table vsearch_output_${seqs}/${seqs}_table-dn-99.qza \
  --o-clustered-sequences vsearch_output_${seqs}/${seqs}_rep-seqs-dn-99.qza \
  --p-perc-identity 0.99

#summarize your filtered/ASV table data
qiime tools export --input-path vsearch_output_${seqs}/${seqs}_demux-filter-stats.qza --output-path vsearch_output_${seqs}/${seqs}

qiime feature-table summarize --i-table vsearch_output_${seqs}/${seqs}_table.qza --o-visualization vsearch_output_${seqs}/${seqs}_table_summary.qzv --verbose

done

### Merging denoised data

# ASV table
qiime feature-table merge \
  --i-tables vsearch_output_${RUN1}/${RUN1}_table-dn-99.qza \
  --i-tables vsearch_output_${RUN2}/${RUN2}_table-dn-99.qza \
  --o-merged-table vsearch_output/table.qza

#summarize
qiime feature-table summarize \
  --i-table vsearch_output/table.qza \
  --o-visualization vsearch_output/table.qzv 
  ##--m-sample-metadata-file sample-metadata.tsv


#qiime feature-table summarize \
#  --i-table table.qza \
#  --o-visualization table.qzv \
#  --m-sample-metadata-file sample-metadata.tsv
#qiime feature-table tabulate-seqs \
#  --i-data rep-seqs.qza \
#  --o-visualization rep-seqs.qzv

### Build quick phylogeny with FastTree
#https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2018.8)
#https://docs.qiime2.org/2018.11/tutorials/moving-pictures/

mkdir phylogeny

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

qiime feature-classifier classify-sklearn \
  --i-classifier /homedir/galati/data/metab/ITS/classifier/unite-ver7-dynamic-classifier-01.12.2017.qza \
  --i-reads vsearch_output/representative_sequences.qza \
  --o-classification taxonomy/ITS_taxonomy.qza \
  --p-n-jobs ${NSLOTS} \
  --verbose

qiime metadata tabulate \
  --m-input-file taxonomy/ITS_taxonomy.qza \
  --o-visualization taxonomy/ITS_taxonomy.qzv

# necessite metadata
# qiime taxa barplot \
#  --i-table vsearch_output/table.qza \
#  --i-taxonomy taxonomy/${DB}_taxonomy.qza \
#  --o-visualization taxonomy/${DB}_taxa-bar-plots.qzv \
#  --m-metadata-file metadata.tsv 

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
