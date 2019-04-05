#!/bin/bash

#$ -q bigmem.q
#$ -N q2_deblur_ITS
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 25
#$ -l mem_free=16G
#$ -V
#$ -cwd

module purge
module load system/conda/5.1.0
source activate qiime2-2018.11

# JOB BEGIN

mkdir /homedir/galati/data/ITS_primer_trimmed2_analysis/
mv /homedir/galati/data/ITS_primer_trimmed2/*.txt /homedir/galati/data/ITS_primer_trimmed2_analysis/
mv /homedir/galati/data/ITS_primer_trimmed2/QC/ /homedir/galati/data/ITS_primer_trimmed2_analysis/

rm -r dada2_output
rm -r deblur_output
rm -r phylogeny
rm -r taxonomy

IN=/homedir/galati/data

RUN1=ITS_primer_trimmed2
RUN2=Mock_S280

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

mkdir deblur_output

for seqs in ${RUN1} ${RUN2}
do

### Quality filter

qiime quality-filter q-score-joined \
 --i-demux ${IN}/${seqs}_demux-joined.qza \
 --o-filtered-sequences deblur_output_${seqs}/${seqs}_demux-joined-filtered.qza \
 --o-filter-stats deblur_output_${seqs}/${seqs}_demux-joined-filter-stats.qza

qiime deblur denoise-other \
  --i-demultiplexed-seqs deblur_output_${seqs}/${seqs}_demux-joined-filtered.qza \
  --p-trim-length 220 \
  --o-representative-sequences deblur_output_${seqs}/${seqs}_rep-seqs-deblur.qza \
  --o-table deblur_output_${seqs}/${seqs}_table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur_output_${seqs}/${seqs}_deblur-stats.qza


### Viewing denoising stats
qiime metadata tabulate \
  --m-input-file deblur_output_${seqs}/${seqs}_demux-joined-filter-stats.qza \
  --o-visualization deblur_output_${seqs}/${seqs}_demux-filter-stats.qzv


qiime deblur visualize-stats \
  --i-deblur-stats deblur_output_${seqs}/${seqs}_deblur-stats.qza \
  --o-visualization deblur_output_${seqs}/${seqs}_deblur-joined-stats.qzv

#summarize your filtered/ASV table data
qiime tools export --input-path deblur_output_${seqs}/${seqs}_demux-filter-stats.qza --output-path deblur_output_${seqs}/${seqs}

qiime feature-table summarize --i-table deblur_output_${seqs}/${seqs}_table-deblur.qza --o-visualization deblur_output_${seqs}/${seqs}_table_summary.qzv --verbose

done

### Merging denoised data

# ASV table
qiime feature-table merge \
  --i-tables deblur_output_${RUN1}/${RUN1}_table-deblur.qza \
  --i-tables deblur_output_${RUN2}/${RUN2}_table-deblur.qza \
  --o-merged-table deblur_output/table.qza

# Representative sequences
qiime feature-table merge-seqs \
  --i-data deblur_output_${RUN1}/${RUN1}_representative_sequences.qza \
  --i-data deblur_output_${RUN2}/${RUN2}_representative_sequences.qza \
  --o-merged-data deblur_output/representative_sequences.qza

# Denoising Stats

cat deblur_output_${RUN1}/${RUN1}/stats.tsv deblur_output_${RUN2}/${RUN2}/stats.tsv \
    > deblur_output/stats.tsv

#cannot
#qiime feature-table merge \
#  --i-tables deblur_output/${RUN1}_denoising_stats.qza  \
#  --i-tables deblur_output/${RUN2}_denoising_stats.qza  \
#  --o-merged-table deblur_output/denoising_stats.qza

#summarize
qiime feature-table summarize \
  --i-table deblur_output/table.qza \
  --o-visualization deblur_output/table.qzv 
  ##--m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data deblur_output/representative_sequences.qza\
  --o-visualization deblur_output/representative_sequences.qzv

#cannot
#qiime feature-table summarize --i-table deblur_output/denoising_stats.qza \
#  --o-visualization deblur_output/denoising_stats.qzv

#export
#qiime tools export deblur_output/denoising_stats.qza --output-path deblur_output

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
  --i-sequences deblur_output/representative_sequences.qza \
  --o-alignment phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree phylogeny/unrooted-tree.qza \
  --o-rooted-tree phylogeny/rooted-tree.qza \
  --verbose

#qiime tools export --input-path deblur_output/deblur_table_filt.qza --output-path deblur_output_exported
#qiime tools export --input-path deblur_output/rep_seqs_filt.qza --output-path deblur_output_exported

### Assign Taxonomy
# loop to test various taxonomic database - pour toi laisser juste silva123 - 
# https://www.dropbox.com/s/5tckx2vhrmf3flp/silva-132-99-nb-classifier.qza?dl=0


mkdir taxonomy

qiime feature-classifier classify-sklearn \
  --i-classifier /homedir/galati/data/classifier/unite-ver7-dynamic-classifier-01.12.2017.qza \
  --i-reads deblur_output/representative_sequences.qza \
  --o-classification taxonomy/ITS_taxonomy.qza \
  --p-n-jobs ${NSLOTS} \
  --verbose

qiime metadata tabulate \
  --m-input-file taxonomy/ITS_taxonomy.qza \
  --o-visualization taxonomy/ITS_taxonomy.qzv

# necessite metadata
# qiime taxa barplot \
#  --i-table deblur_output/table.qza \
#  --i-taxonomy taxonomy/${DB}_taxonomy.qza \
#  --o-visualization taxonomy/${DB}_taxa-bar-plots.qzv \
#  --m-metadata-file metadata.tsv 

qiime tools export --input-path taxonomy/ITS_taxonomy.qza --output-path taxonomy
mv taxonomy/taxonomy.tsv taxonomy/ITS_taxonomy.tsv


### Exporting and modifying BIOM tables

#Creating a TSV BIOM table
qiime tools export --input-path deblur_output/table.qza --output-path export
biom convert -i export/feature-table.biom -o export/ASV-table.biom.tsv --to-tsv

cp export/ASV-table.biom.tsv export/feature-table.biom.tsv

#https://askubuntu.com/questions/20414/find-and-replace-text-within-a-file-using-commands
# ASV table for phyloseq

sed -i "s/#OTU ID/OTUID/g" export/ASV-table.biom.tsv
sed -i "1d" export/ASV-table.biom.tsv

sed -i "s/#OTU ID/#OTUID/g" export/feature-table.biom.tsv

#Export Taxonomy
qiime tools export --input-path /homedir/galati/data/classifier/unite-ver7-dynamic-classifier-01.12.2017.qza --output-path export

biom add-metadata -i export/ASV-table.biom.tsv -o export/ASV-table-silva-132-taxonomy.biom \
  --observation-metadata-fp export/unite-ver7-dynamic_taxonomy.tsv \
  --sc-separated taxonomy
biom convert -i export/ASV-table-unite-ver7-dynamic-taxonomy.biom -o export/ASV-table-unite-ver7-dynamic-taxonomy.biom.tsv --to-tsv

#Export ASV seqs
qiime tools export --input-path deblur_output/representative_sequences.qza --output-path export

#Export Tree
qiime tools export --input-path phylogeny/unrooted-tree.qza --output-path export
mv export/tree.nwk export/unrooted-tree.nwk
qiime tools export --input-path phylogeny/rooted-tree.qza --output-path export
mv export/tree.nwk export/rooted-tree.nwk


# For phyloseq you will need :
#ls export/feature-table.biom.tsv export/taxonomy.tsv export/unrooted-tree.nwk export/rooted-tree.nwk

zip export/export.zip export/* deblur_outpu*/*qzv taxonomy/*.qzv


# JOB END
date

exit 0

