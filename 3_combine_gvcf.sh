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
PICARD=/home/adm-loc/Documents/apps/picard/picard.jar

#FILES contains the BAM file for each sample
FILES_BAM=*.bam

#FILES contains the indexed BAM file for each sample
FILES_BAI=*.bai

#FILES contains the fixed, indexed BAM file for each sample
FILES_FIX=*.fix

#Location of the Reference genome
REF_GEN=/home/adm-loc/Documents/genetic_data/Pfalciparum/fasta/Pfalciparum.genome.fasta

#Location of genetic crosses
GEN_CROSS=/home/adm-loc/Documents/genetic_data/Pfalciparum/known_sites

#Location of VCF files
FILES_VCF=*.vcf

###############################################################################################################
#Create dictionnaries from reference genome with the name of each sample
for f in `ls $FILES_VCF`
do
    vcftools --vcf $f --min-meanDP 5 --recode --recode-INFO-all --out $f.fix
    echo "dictionary created $f"
    rm $f
done

for f in `ls $FILES_VCF`
do 
    mv -- "$f" "${f%.fix.recode.vcf}"
    echo "fixing VCF PROCESSED $f"
done

#Combine VCF files
#Add the vcf for each sample
gatk CombineGVCFs \
    -R $REF_GEN \
    --variant IPC5188_sorted.bam.g.vcf \
    --variant K102_sorted.bam.g.vcf \
    --variant K12_sorted.bam.g.vcf \
    --variant K1_sorted.bam.g.vcf \
    --variant K4_sorted.bam.g.vcf \
    --variant O141A_sorted.bam.g.vcf \
    --variant PA0034-C_sorted.bam.g.vcf \
    --variant PA0040-C_sorted.bam.g.vcf \
    --variant PA0061-C_sorted.bam.g.vcf \
    --variant PA0138-C_sorted.bam.g.vcf \
    --variant PA0144-C_sorted.bam.g.vcf \
    --variant PE0090-C_sorted.bam.g.vcf \
    --variant PE0091-C_sorted.bam.g.vcf \
    --variant PE0105-C_sorted.bam.g.vcf \
    --variant PE0109-C_sorted.bam.g.vcf \
    --variant PE0117-C_sorted.bam.g.vcf \
    --variant PH0210-C_sorted.bam.g.vcf \
    --variant PH0283-C_sorted.bam.g.vcf \
    --variant PH0306-C_sorted.bam.g.vcf \
    --variant PH0364-C_sorted.bam.g.vcf \
    --variant PH0371-C_sorted.bam.g.vcf \
    --variant PV0254-C_sorted.bam.g.vcf \
    --variant PV0256-C_sorted.bam.g.vcf \
    --variant PV0257-C_sorted.bam.g.vcf \
    --variant PV0258-C_sorted.bam.g.vcf \
    --variant PV0259-C_sorted.bam.g.vcf \
    --variant PW0001-C_sorted.bam.g.vcf \
    --variant PW0002-C_sorted.bam.g.vcf \
    --variant PW0003-C_sorted.bam.g.vcf \
    --variant PW0004-C_sorted.bam.g.vcf \
    --variant PW0005-C_sorted.bam.g.vcf \
    --variant PW0006-C_sorted.bam.g.vcf \
    --variant PW0007-C_sorted.bam.g.vcf \
    --variant PW0008-C_sorted.bam.g.vcf \
    --variant PW0009-C_sorted.bam.g.vcf \
    --variant PW0010-CW_sorted.bam.g.vcf \
    --variant S10_sorted.bam.g.vcf \
    --variant S691_sorted.bam.g.vcf \
    --variant U236_sorted.bam.g.vcf \
    -O combine.vcf
echo "CombineGVCFs PROCESSED"

#Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
gatk GenotypeGVCFs \
    -R $REF_GEN \
    -V combine.vcf \
    --max-alternate-alleles 6 \
    -O calling_GVCF.vcf
echo "GenotypeGVCFs PROCESSED"