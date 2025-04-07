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
    $(for f in *.g.vcf; do echo "--variant $f"; done) \
    -O combine.vcf
echo "CombineGVCFs PROCESSED"


#Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
gatk GenotypeGVCFs \
    -R $REF_GEN \
    -V combine.vcf \
    --max-alternate-alleles 6 \
    -O calling_GVCF.vcf
echo "GenotypeGVCFs PROCESSED"