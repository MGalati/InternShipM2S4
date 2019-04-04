WORK IN PROGRESS

#!/bin/bash

#$ -q normal.q
#$ -N q2_vsearch
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 16
#$ -l mem_free=16G
#$ -V
#$ -cwd

module purge
module load system/conda/5.1.0
source activate qiime2-2018.11

# JOB BEGIN

mkdir /homedir/galati/data/16S_primer_trimmed2_analysis/
mv /homedir/galati/data/16S_primer_trimmed2/*.txt /homedir/galati/data/16S_primer_trimmed2_analysis/
mv /homedir/galati/data/16S_primer_trimmed2/QC/ /homedir/galati/data/16S_primer_trimmed2_analysis/
mv /homedir/galati/mock/analysis/16S/pair/Mock_S280/qc/ /homedir/galati/mock/analysis/16S/pair/Mock_S280_qc/
mv /homedir/galati/mock/analysis/16S/pair/Mock_S280 /homedir/galati/data/

rm -r dada2_output
rm -r phylogeny
rm -r taxonomy

IN=/homedir/galati/data

RUN1=16S_primer_trimmed2
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

for seqs in ${RUN1} ${RUN2}
do

### Quality filter

qiime quality-filter q-score \
 --i-demux ${IN}/${seqs}_demux.qza \
 --o-filtered-sequences dada2_output_${seqs}/${seqs}_demux-filtered.qza \
 --o-filter-stats dada2_output_${seqs}/${seqs}_demux-filter-stats.qza


qiime deblur denoise-16S \
  --i-demultiplexed-seqs dada2_output_${seqs}/${seqs}_demux-filtered.qza \
  --p-trim-length 220 \
  --o-representative-sequences dada2_output_${seqs}/${seqs}_rep-seqs-deblur.qza \
  --o-table dada2_output_${seqs}/${seqs}_table-deblur.qza \
  --p-sample-stats \
  --o-stats dada2_output_${seqs}/${seqs}_deblur-stats.qza


### Viewing denoising stats
qiime metadata tabulate \
  --m-input-file dada2_output_${seqs}/${seqs}_demux-filter-stats.qza \
  --o-visualization dada2_output_${seqs}/${seqs}_demux-filter-stats.qzv


qiime deblur visualize-stats \
  --i-deblur-stats dada2_output_${seqs}/${seqs}_deblur-stats.qza \
  --o-visualization dada2_output_${seqs}/${seqs}_deblur-stats.qzv

#summarize your filtered/ASV table data
qiime tools export --input-path dada2_output_${seqs}/${seqs}_demux-filter-stats.qza --output-path dada2_output_${seqs}/${seqs}

qiime feature-table summarize --i-table dada2_output_${seqs}/${seqs}_table-deblur.qza --o-visualization dada2_output_${seqs}/${seqs}_table_summary.qzv --verbose

done

### Merging denoised data

mkdir dada2_output
# ASV table
qiime feature-table merge \
  --i-tables dada2_output_${RUN1}/${RUN1}_table-deblur.qza \
  --i-tables dada2_output_${RUN2}/${RUN2}_table-deblur.qza \
  --o-merged-table dada2_output/table.qza

# Representative sequences
qiime feature-table merge-seqs \
  --i-data dada2_output_${RUN1}/${RUN1}_representative_sequences.qza \
  --i-data dada2_output_${RUN2}/${RUN2}_representative_sequences.qza \
  --o-merged-data dada2_output/representative_sequences.qza

# Denoising Stats

cat dada2_output_${RUN1}/${RUN1}/stats.tsv dada2_output_${RUN2}/${RUN2}/stats.tsv \
    > dada2_output/stats.tsv

#cannot
#qiime feature-table merge \
#  --i-tables dada2_output/${RUN1}_denoising_stats.qza  \
#  --i-tables dada2_output/${RUN2}_denoising_stats.qza  \
#  --o-merged-table dada2_output/denoising_stats.qza

#summarize
qiime feature-table summarize \
  --i-table dada2_output/table.qza \
  --o-visualization dada2_output/table.qzv 
  ##--m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data dada2_output/representative_sequences.qza\
  --o-visualization dada2_output/representative_sequences.qzv

#cannot
#qiime feature-table summarize --i-table dada2_output/denoising_stats.qza \
#  --o-visualization dada2_output/denoising_stats.qzv

#export
#qiime tools export dada2_output/denoising_stats.qza --output-path dada2_output

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
  --i-sequences dada2_output/representative_sequences.qza \
  --o-alignment phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree phylogeny/unrooted-tree.qza \
  --o-rooted-tree phylogeny/rooted-tree.qza \
  --verbose

#qiime tools export --input-path dada2_output/deblur_table_filt.qza --output-path dada2_output_exported
#qiime tools export --input-path dada2_output/rep_seqs_filt.qza --output-path dada2_output_exported

### Assign Taxonomy
# loop to test various taxonomic database - pour toi laisser juste silva123 - 
# https://www.dropbox.com/s/5tckx2vhrmf3flp/silva-132-99-nb-classifier.qza?dl=0


mkdir taxonomy

qiime feature-classifier classify-sklearn \
  --i-classifier /homedir/galati/data/classifier/silva-132-99-nb-classifier.qza \
  --i-reads dada2_output/representative_sequences.qza \
  --o-classification taxonomy/16S_taxonomy.qza \
  --p-n-jobs ${NSLOTS} \
  --verbose

qiime metadata tabulate \
  --m-input-file taxonomy/16S_taxonomy.qza \
  --o-visualization taxonomy/16S_taxonomy.qzv

# necessite metadata
# qiime taxa barplot \
#  --i-table dada2_output/table.qza \
#  --i-taxonomy taxonomy/${DB}_taxonomy.qza \
#  --o-visualization taxonomy/${DB}_taxa-bar-plots.qzv \
#  --m-metadata-file metadata.tsv 

qiime tools export --input-path taxonomy/16S_taxonomy.qza --output-path taxonomy
mv taxonomy/taxonomy.tsv taxonomy/16S_taxonomy.tsv


### Exporting and modifying BIOM tables

#Creating a TSV BIOM table
qiime tools export --input-path dada2_output/table.qza --output-path export
biom convert -i export/feature-table.biom -o export/ASV-table.biom.tsv --to-tsv

cp export/ASV-table.biom.tsv export/feature-table.biom.tsv

#https://askubuntu.com/questions/20414/find-and-replace-text-within-a-file-using-commands
# ASV table for phyloseq

sed -i "s/#OTU ID/OTUID/g" export/ASV-table.biom.tsv
sed -i "1d" export/ASV-table.biom.tsv

sed -i "s/#OTU ID/#OTUID/g" export/feature-table.biom.tsv

#Export Taxonomy
# idem taxonomy
# loop to test various taxonomic database - pour toi laisser juste silva123 - 
# https://www.dropbox.com/s/5tckx2vhrmf3flp/silva-132-99-nb-classifier.qza?dl=0

for DB in silva_123
do

qiime tools export --input-path taxonomy/${DB}_taxonomy.qza --output-path export
mv export/taxonomy.tsv export/${DB}_taxonomy.tsv
#sed -i "s/Feature ID/OTUID/g" export/${DB}_taxonomy.tsv
sed -i "1s/.*/#OTUID\ttaxonomy\tconfidence/" export/${DB}_taxonomy.tsv

##Add taxonomy to biom table
##https://forum.qiime2.org/t/exporting-and-modifying-biom-tables-e-g-adding-taxonomy-annotations/3630

#cp export/taxonomy.tsv export/biom-taxonomy.tsv
##Change the first line of biom-taxonomy.tsv (i.e. the header) to this:
##OTUID  taxonomy  confidence
#sed -i "1s/.*/#OTUID\ttaxonomy\tconfidence/" export/biom-taxonomy.tsv

biom add-metadata -i export/ASV-table.biom.tsv  -o export/ASV-table-${DB}-taxonomy.biom \
  --observation-metadata-fp export/${DB}_taxonomy.tsv \
  --sc-separated taxonomy
biom convert -i export/ASV-table-${DB}-taxonomy.biom -o export/ASV-table-${DB}-taxonomy.biom.tsv --to-tsv
done

#Export ASV seqs
qiime tools export --input-path dada2_output/representative_sequences.qza --output-path export

#Export Tree
qiime tools export --input-path phylogeny/unrooted-tree.qza --output-path export
mv export/tree.nwk export/unrooted-tree.nwk
qiime tools export --input-path phylogeny/rooted-tree.qza --output-path export
mv export/tree.nwk export/rooted-tree.nwk


# For phyloseq you will need :
#ls export/feature-table.biom.tsv export/taxonomy.tsv export/unrooted-tree.nwk export/rooted-tree.nwk

zip export/export.zip export/* dada2_outpu*/*qzv taxonomy/*.qzv

#locally
#scp florentin2@172.28.30.116:/media/DataDrive05/Flo/ETH/export/export.zip .
#
#otu <- read.table(Sile = "feature-table.biom.tsv", header = TRUE)
#tax <- read.table(Sile = "taxonomy.tsv", sep = '\t', header = TRUE)
#merged_Sile <- merge(otu, tax, by.x = c("OTUID"), by.y=c("OTUID"))

#OTU=otu_table(as.matrix(read.csv("otu_matrix.csv", sep=",", row.names=1)), taxa_are_rows = TRUE)
#TAX=tax_table(as.matrix(read.csv("taxonomy.csv", sep=",", row.names=1)))
#TREE = read_tree("rooted-tree.nwk")
#Unifrac requires a rooted tree for calculation. FastUnifrac online selects an arbitrary root when an unrooted tree is uploaded.
#https://github.com/joey711/phyloseq/issues/235

#physeq = phyloseq(OTU, TAX, META, TREE)


# JOB END
date

exit 0
