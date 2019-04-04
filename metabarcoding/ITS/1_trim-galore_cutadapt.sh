#!/bin/bash

#$ -q short.q
#$ -N trim-galore_cutadapt_ITS
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 12
#$ -l mem_free=12G
#$ -V
#$ -cwd
#$ -V

module purge
module load system/conda/5.1.0

# conda create --name trimgalore
# conda install -c bioconda trim-galore 
source activate trimgalore

echo "Récupération des noms d'échantillons"

rm sample_names.tsv
rm /homedir/galati/data/ITS/Undetermined*

MOCK=/homedir/galati/mock/analysis/ITS/pair/ITS_mock24
RAW_ITS=/homedir/galati/data/ITS

#for f in ${MOCK}
#do
#ls ${f}/*R1*gz | cut -d/ -f8 |cut -d_ -f1,2 >> sample_names.tsv
#done

for f in ${RAW_ITS}
do
ls ${f}/*R1*gz | cut -d/ -f6 |cut -d_ -f1,2 >> sample_names.tsv
done

echo "Suppression des adapters ITS"
TRIM_ITS=/homedir/galati/data/ITS_adapter_test
mkdir ${TRIM_ITS}

for NAME in `awk '{print $1}' sample_names.tsv`
do
    echo "Trim galore Sample $NAME Lancement !"
trim_galore --paired -q 0 --nextera --length 0 ${RAW_ITS}/${NAME}_L001_R1_001.fastq.gz ${RAW_ITS}/${NAME}_L001_R2_001.fastq.gz -o ${TRIM_ITS}
    echo "Trim galore Sample ${prefix_ITS} Fin !"
    
done

echo "Renommage des fichiers au format Casava et à l'extension .fastq.gz"
ls ${TRIM_ITS}
#for f in ${TRIM_ITS}/*.fq.gz ; do mv $f $(echo $f | cut -d "_" -f1-5).fastq.gz ; done ;
rename _val_1.fq.gz .fastq.gz ${TRIM_ITS}/*gz
rename _val_2.fq.gz .fastq.gz ${TRIM_ITS}/*gz
ls ${TRIM_ITS}

echo "Suppresison des fichiers de rapport"
# rm ${TRIM_ITS}/*.txt

IN=/homedir/galati/data/ITS_adapter_test/
OUT=/homedir/galati/data/ITS_primer_trimmed2/

mkdir ${OUT}

for NAME in `awk '{print $1}' sample_names.tsv`
do
    echo "Suppression de primers"
    ls ${IN}${NAME}_L001_R1_001.fastq.gz ${IN}${NAME}_L001_R2_001.fastq.gz
cutadapt \
    --pair-filter any \
    --no-indels \
    --discard-untrimmed \
    -g NNNNGTGAATCATCGAATCTTTGAA \
    -G NNNNTCCTCCGCTTATTGATATGC \
    -o ${OUT}${NAME}"_L001_R1_001.fastq.gz" \
    -p ${OUT}${NAME}"_L001_R2_001.fastq.gz" \
    ${IN}${NAME}*R1*.fastq.gz ${IN}${NAME}*R2*.fastq.gz \
	> ${OUT}${NAME}"_cutadapt_log.txt"
	
#    --max-n 0 \
#    --minimum-length 175 \ 
# Dada2 peut le faire plus tard mais pas sous Qiime2 ...

echo "Fin de la suppression de primers"
done

source deactivate trimgalore
source activate fastqc_multiqc

mkdir ${OUT}'QC/'
fastqc -t ${NSOLTS} ${OUT}/*fastq* -o ${OUT}'QC/'
multiqc ${OUT}'QC/' -o ${OUT}'QC/'


# JOB END
date

exit 0
