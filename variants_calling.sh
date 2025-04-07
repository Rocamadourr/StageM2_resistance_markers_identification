#! /bin/bash

#Script Romain Coppee
#Creation data: 09/02/2020
#Last modification: 09/07/2020

####--------General Goal: Produce VCF files containing only high-quality SNPs from whole genome sequencing data

###############################################################################################################
#####-------Preparation/location of softwares and materials
#Adding samtools/bcftools to the PATH environment
#export PATH="$PATH:/usr/bin/bcftools-1.9"
#export PATH="$PATH:/usr/bin/samtools-1.9"
#export PATH="$PATH:/usr/bin/htslib-1.9"
#source ~/.profile

#Adding GATK to the PATH environment
#export PATH=/home/virologie/Documents/gatk-4.1.8.1/:$PATH

#location of PICARD software
PICARD=c:\Users\Utilisateur\Downloads\stage_noe\    picard.jar

#FILES contains the BAM file for each sample
FILES_BAM=*.bam

#FILES contains the indexed BAM file for each sample
FILES_BAI=*.bai

#FILES contains the fixed, indexed BAM file for each sample
FILES_FIX=*.fix

#Location of the Reference genome
REF_GEN=C:\Users\Utilisateur\Downloads\stage_noe\TEST_NGS\Pfalciparum.genome.fasta

#Location of genetic crosses
GEN_CROSS=/home/adm-loc/Documents/genetic_data/Pfalciparum/known_sites


###############################################################################################################
#####-------Goal 1: Complete Header to fix BAM files, then index them



#Mark Duplicates with PicardTools
for f in `ls $FILES_BAM`
do
    java -jar $PICARD MarkDuplicates \
        -INPUT $f \
        -REMOVE_DUPLICATES true \
        -VALIDATION_STRINGENCY LENIENT \
        -AS true \
        -METRICS_FILE metrics \
        -CREATE_INDEX true \
        -OUTPUT $f.fix
    echo "Mark Duplicates PROCESSED $f"
done

#first remove all .bam files, and rename the fixed .bam.fix to .bam
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing mark duplicates PROCESSED $f"
done



#complete the header of each bam file using the AddOrReplaceReadGroups from Picard 
for f in `ls $FILES_BAM`
do
    file_name="${f##*/}"
    onlyFileName="${file_name%.*}" #remove.bam for SM information
    onlyFileName="${onlyFileName%.*}" #remove .sorted SM information
    java -jar $PICARD AddOrReplaceReadGroups \
        -I $f \
        -O $f.fix \
        -LB WGA \
        -PL illumina \
        -PU NA \
        -SM $onlyFileName
    echo "AddOrReplaceReadGroups PROCESSED $f"
done


#first remove all .bam files, and rename the fixed .bam.fix to .bam
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing AddOrReplaceReadGroups PROCESSED $f"
done

#index each fixed bam file using samtools
for f in `ls $FILES_BAM`
do
    samtools index $f
    echo "index PROCESSED $f"
done


#Recalibrate base qualities in a .bam file so that quality metrics match actual observed error rates
for f in `ls $FILES_BAM`
do
    gatk BaseRecalibrator \
        -R $REF_GEN \
        -I $f \
        --known-sites $GEN_CROSS/3d7_hb3.combined.final.vcf.gz \
        --known-sites $GEN_CROSS/hb3_dd2.combined.final.vcf.gz \
        --known-sites $GEN_CROSS/7g8_gb4.combined.final.vcf.gz \
        -O $f.pass1.table 
    echo "BaseRecalibrator PROCESSED $f"
done

#Apply BQSR to input BAM file
for f in `ls $FILES_BAM`
do
    gatk ApplyBQSR \
        -R $REF_GEN \
        -I $f \
        -bqsr $f.pass1.table \
        -O $f.fix
done

#first remove all .bam files, and rename the fixed .bam.fix to .bam
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing PrintReads PROCESSED $f"
done

#remove the previous indexed bam files and create new indexed files
rm $FILES_BAI
for f in `ls $FILES_BAM`
do
    samtools index $f
    echo "index PROCESSED $f"
done

#generate a complete variant calling from each fixed bam file
for f in `ls $FILES_BAM`
do
    gatk HaplotypeCaller \
        -R $REF_GEN \
        -I $f \
        -ERC GVCF \
        --native-pair-hmm-threads 8 \
        --max-alternate-alleles 6 \
        -O $f.g.vcf
    echo "HaplotypeCaller PROCESSED $f"
done
