#!/bin/bash

#$ -q bigmem.q
#$ -N a-p_16S
#$ -M mathias.galati@cirad.fr
#$ -pe parallel_smp 10
#$ -l mem_free=8G
#$ -V
#$ -cwd

module purge
module load system/conda/5.1.0

# conda create --name trimgalore
# conda install -c bioconda trim-galore 
source activate trimgalore

echo "Récupération des noms d'échantillons"

MOCK=/homedir/galati/data/metab/16S/MOCK/Mock_S280
RAW_16S=/homedir/galati/data/metab/16S/RAW

#for f in ${MOCK}
#do
#ls ${f}/*R1*gz | cut -d/ -f8 |cut -d_ -f1,2 >> sample_names.tsv
#done

for f in ${RAW_16S}
do
ls ${f}/*R1*gz | cut -d/ -f8 |cut -d_ -f1,2 >> sample_names.tsv
done

echo "Suppression des adapters 16S"
TRIM_16S=/homedir/galati/data/metab/16S/TRIM
mkdir ${TRIM_16S}

trim_galore --paired -q 0 --nextera --length 0 ${MOCK}/Mock_S280_L001_R1_001.fastq.gz ${MOCK}/Mock_S280_L001_R2_001.fastq.gz -o ${TRIM_16S}

for NAME in `awk '{print $1}' sample_names.tsv`
do
    echo "Trim galore Sample $NAME Lancement !"
trim_galore --paired -q 0 --nextera --length 0 ${RAW_16S}/${NAME}_L001_R1_001.fastq.gz ${RAW_16S}/${NAME}_L001_R2_001.fastq.gz -o ${TRIM_16S}
    echo "Trim galore Sample ${prefix_16S} Fin !"
    
done

echo "Renommage des fichiers au format Casava et à l'extension .fastq.gz"
ls ${TRIM_16S}
#for f in ${TRIM_16S}/*.fq.gz ; do mv $f $(echo $f | cut -d "_" -f1-5).fastq.gz ; done ;
rename _val_1.fq.gz .fastq.gz ${TRIM_16S}/*gz
rename _val_2.fq.gz .fastq.gz ${TRIM_16S}/*gz
ls ${TRIM_16S}

echo "Suppresison des fichiers de rapport"
# rm ${TRIM_16S}/*.txt

IN=/homedir/galati/data/metab/16S/TRIM/
OUT=/homedir/galati/data/metab/16S/PRIM/

mkdir ${OUT}

for NAME in `awk '{print $1}' sample_names.tsv`
do
    echo "Suppression de primers"
    ls ${IN}${NAME}_L001_R1_001.fastq.gz ${IN}${NAME}_L001_R2_001.fastq.gz
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
	
#    --max-n 0 \
#    --minimum-length 175 \ 
# Dada2 peut le faire plus tard mais pas sous Qiime2 ...

# Il faut que les noms de fichiers respectent le format Casava
# par exemple L2S357_15_L001_R1_001.fastq.gz

echo "Fin de la suppression de primers"
done

source deactivate trimgalore
source activate fastqc_multiqc

mkdir ${OUT}'QC/'
fastqc -t ${NSOLTS} ${OUT}*fastq* -o ${OUT}'QC/'
multiqc ${OUT}'QC/' -o ${OUT}'QC/'

mkdir /homedir/galati/data/metab/16S/PRIM_qc/
mv /homedir/galati/data/metab/16S/PRIM/qc/ ../PRIM_qc
mv /homedir/galati/data/metab/16S/PRIM/*.txt ../PRIM_qc

# JOB END
date

exit 0
