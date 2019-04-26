#!/bin/bash

#$ -q bigmem.q
#$ -N Sdeblur_ITS
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 6
#$ -l mem_free=6G
#$ -V
#$ -cwd

module purge
module load system/conda/5.1.0
source activate qiime2-2018.11

# JOB BEGIN

IN=/homedir/galati/data/metab/ITS

RUN1=SUB
RUN2=ITS_mock24

#mkdir /homedir/galati/data/metab/ITS/deblur
deblur=/homedir/galati/data/metab/ITS/deblur
"""
for seqs in ${RUN1} ${RUN2}
do
echo 'Import'
qiime tools import --type SampleData[JoinedSequencesWithQuality] \
                   --input-path ${IN}/${seqs} \
                   --output-path ${deblur}/${seqs}_reads.qza \
                   --input-format CasavaOneEightSingleLanePerSampleDirFmt 

echo 'Check this artifact to make sure that QIIME now recognizes your data'
qiime tools peek ${deblur}/${seqs}_reads.qza 

echo 'Initial sequence quality control'
qiime demux summarize \
  --i-data ${deblur}/${seqs}_reads.qza  \
  --o-visualization ${deblur}/${seqs}_reads.qzv  \
  --verbose
done

mkdir deblur_output

for seqs in ${RUN1} ${RUN2}
do

echo 'Joining reads'
qiime vsearch join-pairs \
  --i-demultiplexed-seqs ${deblur}/${seqs}_reads.qza \
  --o-joined-sequences ${deblur}/${seqs}_reads-joined.qza

echo 'Summary ${seqs}'
qiime demux summarize \
  --i-data ${deblur}/${seqs}_reads-joined.qza \
  --o-visualization ${deblur}/${seqs}_reads-joined.qzv

mkdir deblur_output_${seqs}
echo 'Quality filter'
qiime quality-filter q-score-joined \
 --i-demux ${deblur}/${seqs}_reads-joined.qza \
 --o-filtered-sequences deblur_output_${seqs}/${seqs}_demux-joined-filtered.qza \
 --o-filter-stats deblur_output_${seqs}/${seqs}_demux-joined-filter-stats.qza

#echo 'Import RefSeq'
#qiime tools import \
#  --input-path /homedir/galati/data/metab/ITS/RefSeq/UNITE_sh_refs_qiime_ver8_dynamic_02.02.2019.fasta \
#  --output-path /homedir/galati/data/metab/ITS/RefSeq/UNITE_sh_refs_qiime_ver8_dynamic_02.02.2019.qza \
#  --type 'FeatureData[Sequence]'

echo 'Denoise deblur'
qiime deblur denoise-other \
  --i-demultiplexed-seqs deblur_output_${seqs}/${seqs}_demux-joined-filtered.qza \
  --i-reference-seqs /homedir/galati/data/metab/ITS/RefSeq/UNITE_sh_refs_qiime_ver8_dynamic_02.02.2019.qza \
  --p-trim-length 220 \
  --o-representative-sequences deblur_output_${seqs}/${seqs}_rep-seqs-deblur.qza \
  --o-table deblur_output_${seqs}/${seqs}_table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur_output_${seqs}/${seqs}_deblur-stats.qza

echo 'Viewing denoising stats'
qiime metadata tabulate \
  --m-input-file deblur_output_${seqs}/${seqs}_demux-joined-filter-stats.qza \
  --o-visualization deblur_output_${seqs}/${seqs}_demux-joined-filter-stats.qzv

echo 'Summarize your filtered/ASV table data'
qiime tools export --input-path deblur_output_${seqs}/${seqs}_demux-joined-filter-stats.qza --output-path deblur_output_${seqs}/${seqs}

echo 'Summarize your deblur data'
qiime tools export --input-path deblur_output_${seqs}/${seqs}_deblur-stats.qza --output-path deblur_output_${seqs}/${seqs}

echo 'Feature table summarize'
qiime feature-table summarize --i-table deblur_output_${seqs}/${seqs}_table-deblur.qza --o-visualization deblur_output_${seqs}/${seqs}_table_summary.qzv --verbose

done
echo 'Fin de la bouche vsearch'

echo 'Merging denoised data'

echo 'ASV table'
qiime feature-table merge \
  --i-tables deblur_output_${RUN1}/${RUN1}_table-deblur.qza \
  --i-tables deblur_output_${RUN2}/${RUN2}_table-deblur.qza \
  --o-merged-table deblur_output/table.qza
"""
echo '******************************************************************************************************'
echo 'Representative sequences'
qiime feature-table merge-seqs \
  --i-data deblur_output_${RUN1}/${RUN1}_rep-seqs-deblur.qza \
  --i-data deblur_output_${RUN2}/${RUN2}_rep-seqs-deblur.qza \
  --o-merged-data deblur_output/representative_sequences.qza

echo 'Denoising Stats'

cat deblur_output_${RUN1}/${RUN1}/stats.tsv deblur_output_${RUN2}/${RUN2}/stats.tsv \
    > deblur_output/stats.tsv

echo 'summarize'
qiime feature-table summarize \
  --i-table deblur_output/table.qza \
  --o-visualization deblur_output/table.qzv 
  --m-sample-metadata-file /homedir/galati/data/metab/ITS/metadata/sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data deblur_output/representative_sequences.qza\
  --o-visualization deblur_output/representative_sequences.qzv

echo 'export'
qiime tools export deblur_output/denoising_stats.qza --output-path deblur_output

echo 'Summarize'
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file /homedir/galati/data/metab/ITS/metadata/sample-metadata.tsv

echo 'Tabulate sequences'
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

mkdir phylogeny

qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads ${NSLOTS} \
  --i-sequences deblur_output/representative_sequences.qza \
  --o-alignment phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree phylogeny/unrooted-tree.qza \
  --o-rooted-tree phylogeny/rooted-tree.qza \
  --verbose

mkdir taxonomy
echo 'clasif'
qiime feature-classifier classify-sklearn \
  --i-classifier /homedir/galati/data/metab/ITS/classifier/unite-ver7-dynamic-classifier-01.12.2017.qza \
  --i-reads deblur_output/representative_sequences.qza \
  --o-classification taxonomy/ITS_taxonomy.qza \
  --p-n-jobs ${NSLOTS} \
  --verbose

qiime metadata tabulate \
  --m-input-file taxonomy/ITS_taxonomy.qza \
  --o-visualization taxonomy/ITS_taxonomy.qzv

echo 'necessite metadata'
qiime taxa barplot \
  --i-table deblur_output/table.qza \
  --i-taxonomy taxonomy/ITS_taxonomy.qza \
  --o-visualization taxonomy/ITS_taxa-bar-plots.qzv \
  --m-metadata-file /homedir/galati/data/metab/ITS/metadata/metadata.tsv 

echo 'export'
qiime tools export --input-path taxonomy/16S_taxonomy.qza --output-path taxonomy
mv taxonomy/taxonomy.tsv taxonomy/16S_taxonomy.tsv


echo 'Exporting and modifying BIOM tables'
echo 'Creating a TSV BIOM table'
qiime tools export --input-path deblur_output/table.qza --output-path export
biom convert -i export/feature-table.biom -o export/ASV-table.biom.tsv --to-tsv

cp export/ASV-table.biom.tsv export/feature-table.biom.tsv

#https://askubuntu.com/questions/20414/find-and-replace-text-within-a-file-using-commands
# ASV table for phyloseq

sed -i "s/#OTU ID/OTUID/g" export/ASV-table.biom.tsv
sed -i "1d" export/ASV-table.biom.tsv

sed -i "s/#OTU ID/#OTUID/g" export/feature-table.biom.tsv

echo 'Export Taxonomy'
qiime tools export --input-path /homedir/galati/data/metab/ITS/classifier/unite-ver7-dynamic-classifier-01.12.2017.qza --output-path export

biom add-metadata -i export/ASV-table.biom.tsv -o export/ASV-table-silva-132-taxonomy.biom \
  --observation-metadata-fp export/unite-ver7-dynamic_taxonomy.tsv \
  --sc-separated taxonomy
biom convert -i export/ASV-table-unite-ver7-dynamic-taxonomy.biom -o export/ASV-table-unite-ver7-dynamic-taxonomy.biom.tsv --to-tsv

echo 'Export ASV seqs'
qiime tools export --input-path deblur_output/representative_sequences.qza --output-path export

echo 'Export Tree'
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
