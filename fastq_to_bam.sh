#!/bin/bash

# Folder containing FASTQ
folder="C:/Users/Utilisateur/Downloads/stage_noe"

# Reference fasta (MalariaGen)
reference="Pfalciparum.genome.fasta"

# Verification of the reference file
if [ ! -f "$reference" ]; then
    echo "ERROR : '$reference' not found."
    exit 1
fi

# Checking if R1 files exist
r1_count=$(ls ${folder}/*_R1.fastq 2>/dev/null | wc -l)
if [ $r1_count -eq 0 ]; then
    echo "ERROR : no R1 file (*_R1.fastq)found in the folder '$folder'."
    exit 1
fi

# Repertoire de sortie pour les fichiers SAM, BAM et indexes
mkdir -p ${folder}/output

# Loop througt all the FASTQ in the folder
for r1 in ${folder}/*_R1.fastq; do
    r2=${r1/_R1.fastq/_R2.fastq}  # Checking for the corresponding R2 file
    
    if [ -f "$r2" ]; then
        echo "Treating files : $r1 and $r2"

        # Basename for the output files
        base=$(basename $r1 _R1.fastq)
        sam_output="${folder}/output/${base}.sam"
        bam_output="${folder}/output/${base}.bam"
        sorted_bam_output="${folder}/output/${base}_sorted.bam"

        bwa mem -t 12 $reference $r1 $r2 > $sam_output
        samtools view -b $sam_output > $bam_output
        samtools sort $bam_output -o $sorted_bam_output
        samtools index $sorted_bam_output

        echo "Files generated : ${sorted_bam_output} and ${sorted_bam_output}.bai"
        rm $sam_output $bam_output
    else
        echo "ERROR : corresponding file R2 not found for $r1"
    fi
done
