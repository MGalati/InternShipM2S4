#!/bin/bash

#$ -q normal.q
#$ -N trim-galore_cutadapt
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 6
#$ -l mem_free=6G
#$ -V
#$ -cwd
#$ -V


module purge
module load system/conda/5.1.0

# conda create --name trimgalore
# conda install -c bioconda trim-galore 
conda activate trimgalore

echo 'Récupération des noms d'échantillons'

MOCK=/homedir/galati/mock/analysis/16S/pair/Mock_S280
RAW_16S=/homedir/galati/data/16S

for f in ${MOCK} ${RAW_16S}
do
ls ${f}/*R1*gz | cut -d/ -f2 |cut -d_ -f1,2 >> sample_names.tsv
done

echo "Suppression des adapters 16S"
TRIM_16S=/homedir/galati/data/16S_adapter_test/

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

IN=/homedir/galati/data/16S_adapter_test/
OUT=/homedir/galati/data/16S_primer_trimmed/

mkdir $OUT

for NAME in `awk '{print $1}' sample_names.tsv`
do
    echo 'Suppression de primers'
    ls $IN$NAME*R1*gz $IN$NAME*R2*gz
cutadapt \
    --pair-filter any \
    --no-indels \
    --discard-untrimmed \
    -g NNNNCCTACGGGNBGCASCAG \
    -G NNNNGACTACNVGGGTATCTAATCC \
    -o ${OUT}${NAME}"_L001_R1_001.fastq.gz" \
    -p ${OUT}${NAME}"_L001_R2_001.fastq.gz" \
    ${IN}${NAME}*R1*.fastq.gz ${IN}${NAME}*R2*.fastq.gz \
	> ${OUT}${NAME}"_cutadapt_log.txt"
	
#   --max-n 0 \
#   --minimum-length 175 \ 
# Dada2 can do it later but not within qiime2... so we might do it at this step

# rename samples to get RUN info and respect qiime2 input format
# e.g., L2S357_15_L001_R1_001.fastq.gz

echo 'Fin de la suppression de primers'
done

/homedir/constancias/tools/parse_cutadapt_logs.py \
       -i ${OUT}'/'*log.txt

source deactivate trimgalore
source activate fastqc_multiqc

mkdir ${OUT}'QC/'
fastqc -t ${NSOLTS} ${OUT}*fastq* -o ${OUT}'QC/'
multiqc ${OUT}'QC/' -o ${OUT}'QC/'



# JOB END
date

exit 0
