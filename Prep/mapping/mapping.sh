#!/bin/sh

# --- 以下の変数は適宜変更する ---
# mapping.shではutilは使わない
HISAT2="hisat2"
HISAT2BUILD="${HISAT2}-build"
SAMTOOLS="samtools"
STRINGTIE="stringtie"
# --- ---

if test $# -ne 5; then
    echo '1: genome.fasta <STRING>'
    echo '2: RNA_seq_read1.fastq <STRING>'
    echo '3: RNA_seq_read2.fastq <STRING>'
    echo '4: prefix in output file name <STRING>'
    echo '5: the number or threads'
else
    mkdir ${4}_tmp_index
    mkdir /scratch/tmp_mapping_${4}
    # mapping
    ${HISAT2BUILD} -p ${5} ${1} ${4}_tmp_index/${4}_index
    ${HISAT2} -p ${5} -x ${4}_tmp_index/${4}_index -1 ${2} -2 ${3} -S /scratch/tmp_mapping_${4}/${4}.sam --dta --no-discordant --no-mixed
    ${SAMTOOLS} view -bS -@ ${5} /scratch/tmp_mapping_${4}/${4}.sam | ${SAMTOOLS} sort -@ ${5} -o /scratch/tmp_mapping_${4}/${4}.sorted.bam
    
    # gene prediction
    ${STRINGTIE} -p ${5} -o ${4}.gtf /scratch/tmp_mapping_${4}/${4}.sorted.bam
    
    rm -rf ${4}_tmp_index
    rm -rf /scratch/tmp_mapping_${4}
fi
