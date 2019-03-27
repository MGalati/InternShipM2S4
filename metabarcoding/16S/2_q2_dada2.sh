#!/bin/bash

#$ -q long.q
#$ -N q2_dada2
#$ -M florentin.constancias@cirad.fr
#$ -pe parallel_smp 25
#$ -l mem_free=6G
#$ -V
#$ -cwd

module purge
module load system/conda/5.1.0
source activate qiime2-2018.11

# JOB BEGIN

IN=primer_trimmed_fastqs

# remove mv non .fastq.gz files
# rm ${IN}/*txt # we summarized those before
# mv ${IN}/qc ../


# mv fastq.gz
for RUN in BP CP DD IB II LU MC
do
mkdir -p ${IN}/${RUN}
mv ${IN}/${RUN}*gz ${IN}/${RUN}
done

### Import sequences 
# only sequences and properly formatted in the INPUT directory
# in the following directories

RUN1=BP
RUN2=CP 
RUN3=DD
RUN4=IB
RUN5=II
RUN6=LU
RUN7=MC

for seqs in ${RUN1} ${RUN2} ${RUN3} ${RUN4} ${RUN5} ${RUN6} ${RUN7}
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

#qiime tools view ${seqs}_reads.qzv

### Adapter removal
#https://docs.qiime2.org/2018.11/plugins/available/cutadapt/trim-paired/
#https://forum.qiime2.org/t/qiime-cutadapt-trim-paired-failing-silently-in-qiime-2018-6/7468

##qiime cutadapt trim-paired \
##        --i-demultiplexed-sequences ${OUT}/${IN}.qza \
##        --p-cores ${NSLOTS} \
##        --p-front-f "${fwdcut[$i]}" \
##        --p-front-r "${revcut[$i]}" \
##        --p-error-rate 0.1 \
##        --o-trimmed-sequences "$trpth"/"${otpth[$i]}" \
##        --verbose | tee "$trpth"/"${log[$i]}"

### DADA2 workflow
#https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2018.8)
#For multiple run https://docs.qiime2.org/2018.11/tutorials/fmt/

# script run_dada_paired.R, modified to include
# minLen = 175, maxN = 0, 
# dada(derepRs, selfConsist=TRUE, ..., MAX_CONSIST=20)
# vim /home/florentin2/anaconda3/envs/qiime2-2018.11/bin/run_dada_paired.R to include

for seqs in ${RUN1} ${RUN2} ${RUN3} ${RUN4} ${RUN5} ${RUN6} ${RUN7}
do

truncF=0
truncL=220
trimF=0
trimL=0
maxee=2
truncq=10
nreadslearn=10000000
# pour tester

#6000000
#0 - Use all input reads 
#1000000 default
#100000000 working
#1000000000 working MSQ4 & MSQ6
#100000000000 NOT working MSQ4
chim=consensus
#consensu = faster
#pooled

mkdir dada2_output_${seqs}

qiime dada2 denoise-paired --i-demultiplexed-seqs ${IN}/${seqs}_reads.qza \
                           --p-trunc-len-f ${truncF} \
                           --p-trunc-len-r ${truncL} \
                           --p-trim-left-f ${trimF} \
                           --p-trim-left-r ${trimL} \
                           --p-max-ee ${maxee} \
                           --p-trunc-q ${truncq} \
                           --p-n-reads-learn ${nreadslearn} \
                           --p-n-threads ${NSLOTS} \
                           --p-chimera-method ${chim} \
                           --o-representative-sequences dada2_output_${seqs}/${seqs}_representative_sequences.qza \
                           --o-table dada2_output_${seqs}/${seqs}_table.qza \
                           --o-denoising-stats dada2_output_${seqs}/${seqs}_denoising_stats.qza \
                           --verbose
#                          --output-dir dada2_output \

### Viewing denoising stats
qiime metadata tabulate \
  --m-input-file dada2_output_${seqs}/${seqs}_denoising_stats.qza \
  --o-visualization dada2_output_${seqs}/${seqs}_denoising_stats.qza

#summarize your filtered/ASV table data
qiime tools export --input-path dada2_output_${seqs}/${seqs}_denoising_stats.qza --output-path dada2_output_${seqs}/${seqs}

qiime feature-table summarize --i-table dada2_output_${seqs}/${seqs}_table.qza --o-visualization dada2_output_${seqs}/${seqs}_table_summary.qzv --verbose

done

### Merging denoised data

mkdir dada2_output
# ASV table
qiime feature-table merge \
  --i-tables dada2_output_${RUN1}/${RUN1}_table.qza \
  --i-tables dada2_output_${RUN2}/${RUN2}_table.qza \
  --i-tables dada2_output_${RUN3}/${RUN3}_table.qza \
  --i-tables dada2_output_${RUN4}/${RUN4}_table.qza \
  --i-tables dada2_output_${RUN5}/${RUN5}_table.qza \
  --i-tables dada2_output_${RUN6}/${RUN6}_table.qza \
  --i-tables dada2_output_${RUN7}/${RUN7}_table.qza \
  --o-merged-table dada2_output/table.qza

# Representative sequences
qiime feature-table merge-seqs \
  --i-data dada2_output_${RUN1}/${RUN1}_representative_sequences.qza \
  --i-data dada2_output_${RUN2}/${RUN2}_representative_sequences.qza \
  --i-data dada2_output_${RUN3}/${RUN3}_representative_sequences.qza \
  --i-data dada2_output_${RUN4}/${RUN4}_representative_sequences.qza \
  --i-data dada2_output_${RUN5}/${RUN5}_representative_sequences.qza \
  --i-data dada2_output_${RUN6}/${RUN6}_representative_sequences.qza \
  --i-data dada2_output_${RUN7}/${RUN7}_representative_sequences.qza \
  --o-merged-data dada2_output/representative_sequences.qza

# Denoising Stats

cat dada2_output_${RUN1}/${RUN1}/stats.tsv dada2_output_${RUN2}/${RUN2}/stats.tsv \
    dada2_output_${RUN3}/${RUN3}/stats.tsv dada2_output_${RUN4}/${RUN4}/stats.tsv \
    dada2_output_${RUN5}/${RUN5}/stats.tsv dada2_output_${RUN6}/${RUN6}/stats.tsv \
    dada2_output_${RUN7}/${RUN7}/stats.tsv > dada2_output/stats.tsv

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

for DB in silva_132_99_16S_majority_taxonomy_CCTACGGGNBGCASCAG_GACTACNVGGGTATCTAATCC_id0.85_ml350_Ml500 gg-13-8-99-nb-classifier
#silva_132_99_16S_consensus_taxonomy_GTGCCAGCMGCCGCGGTAA_GGACTACHVGGGTWTCTAAT_id0.85_ml248_Ml257 silva_132_99_16S_consensus_taxonomy_GTGCCAGCMGCCGCGGTAA_GGACTACHVGGGTWTCTAAT_id0.8_ml210_Ml310 \
#gg-13-8-99-515-806-nb-classifier gg-13-8-99-nb-classifier \
#silva-132-99-515-806-nb-classifier silva-132-99-nb-classifier \
#silva_132_99_16S_majority_taxonomy_GTGCCAGCMGCCGCGGTAA_GGACTACHVGGGTWTCTAAT_id0.85_ml248_Ml257

do

qiime feature-classifier classify-sklearn \
  --i-classifier /homedir/constancias/db/qiime2/${DB}.qza \
  --i-reads dada2_output/representative_sequences.qza \
  --o-classification taxonomy/${DB}_taxonomy.qza \
  --p-n-jobs ${NSLOTS} \
  --verbose

qiime metadata tabulate \
  --m-input-file taxonomy/${DB}_taxonomy.qza \
  --o-visualization taxonomy/${DB}_taxonomy.qzv

# necessite metadata
# qiime taxa barplot \
#  --i-table dada2_output/table.qza \
#  --i-taxonomy taxonomy/${DB}_taxonomy.qza \
#  --o-visualization taxonomy/${DB}_taxa-bar-plots.qzv \
#  --m-metadata-file metadata.tsv 

qiime tools export --input-path taxonomy/${DB}_taxonomy.qza --output-path taxonomy
mv taxonomy/taxonomy.tsv taxonomy/${DB}_taxonomy.tsv

done

### Exporting and modifying BIOM tables

#Creating a TSV BIOM table
qiime tools export --input-path dada2_output/table.qza --output-path export
biom convert -i export/feature-table.biom -o export/ASV-table.biom.tsv --to-tsv

cp export/ASV-table.biom.tsv export/feature-table.biom.tsv

#https://askubuntu.com/questions/20414/find-and-replace-text-within-a-file-using-commands
# ASV table for phyloseq

sed -i 's/#OTU ID/OTUID/g' export/ASV-table.biom.tsv
sed -i '1d' export/ASV-table.biom.tsv

sed -i 's/#OTU ID/#OTUID/g' export/feature-table.biom.tsv

#Export Taxonomy
# idem taxonomy
# loop to test various taxonomic database - pour toi laisser juste silva123 - 
# https://www.dropbox.com/s/5tckx2vhrmf3flp/silva-132-99-nb-classifier.qza?dl=0

for DB in silva_132_99_16S_majority_taxonomy_CCTACGGGNBGCASCAG_GACTACNVGGGTATCTAATCC_id0.85_ml350_Ml500 gg-13-8-99-nb-classifier
#silva_132_99_16S_consensus_taxonomy_GTGCCAGCMGCCGCGGTAA_GGACTACHVGGGTWTCTAAT_id0.85_ml248_Ml257 silva_132_99_16S_consensus_taxonomy_GTGCCAGCMGCCGCGGTAA_GGACTACHVGGGTWTCTAAT_id0.8_ml210_Ml310 \
#gg-13-8-99-515-806-nb-classifier gg-13-8-99-nb-classifier \
#silva-132-99-515-806-nb-classifier silva-132-99-nb-classifier \
#silva_132_99_16S_majority_taxonomy_GTGCCAGCMGCCGCGGTAA_GGACTACHVGGGTWTCTAAT_id0.85_ml248_Ml257
do

qiime tools export --input-path taxonomy/${DB}_taxonomy.qza --output-path export
mv export/taxonomy.tsv export/${DB}_taxonomy.tsv
#sed -i 's/Feature ID/OTUID/g' export/${DB}_taxonomy.tsv
sed -i '1s/.*/#OTUID\ttaxonomy\tconfidence/' export/${DB}_taxonomy.tsv

##Add taxonomy to biom table
##https://forum.qiime2.org/t/exporting-and-modifying-biom-tables-e-g-adding-taxonomy-annotations/3630

#cp export/taxonomy.tsv export/biom-taxonomy.tsv
##Change the first line of biom-taxonomy.tsv (i.e. the header) to this:
##OTUID  taxonomy  confidence
#sed -i '1s/.*/#OTUID\ttaxonomy\tconfidence/' export/biom-taxonomy.tsv

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


