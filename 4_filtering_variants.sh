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

#location of SNPEFF software
SNPEFF=/home/adm-loc/Documents/apps/snpEff/snpEff.jar

#Location of the Reference genome
REF_GEN=/home/adm-loc/Documents/genetic_data/Pfalciparum/fasta/Pfalciparum.genome.fasta

#FILES contains the BAM file for each sample
FILES_BAM=*.bam

#FILES contains the indexed BAM file for each sample
FILES_BAI=*.bai

#FILES contains the fixed, indexed BAM file for each sample
FILES_FIX=*.fix

#Location of genetic crosses
GEN_CROSS=/home/adm-loc/Documents/genetic_data/Pfalciparum/known_sites

#Location of VCF files
FILES_VCF=*.vcf

#Location of annotation files
FILES_ANNOT=/home/adm-loc/Documents/genetic_data/Pfalciparum/annotations


###############################################################################################################


#####-------Goal 3: Filtering variants (only SNPs)

#Build a recalibration model to score variant quality for filtering purposes
#maxGaussians to change if needed (especially for one or two samples)
gatk VariantRecalibrator \
    -R $REF_GEN \
    -V calling_GVCF.vcf \
    --resource:cross1,known=false,training=true,truth=true,prior=15.0 $GEN_CROSS/3d7_hb3.combined.final.vcf.gz \
    --resource:cross2,known=false,training=true,truth=true,prior=15.0 $GEN_CROSS/hb3_dd2.combined.final.vcf.gz \
    --resource:cross3,known=false,training=true,truth=true,prior=15.0 $GEN_CROSS/7g8_gb4.combined.final.vcf.gz \
    -mode SNP \
    -an QD \
    -an FS \
    -an SOR \
    -an DP \
    --max-gaussians 8 \
    -mq-cap 70 \
    -O recal_GVCF.vcf \
    --tranches-file recal_GVCF.tranches \
    --rscript-file recal_GVCF.plots.R

echo "VariantRecalibrator PROCESSED"

#Apply a score cutoff to filter variants based on a recalibration table
gatk ApplyVQSR \
    -R $REF_GEN \
    -V calling_GVCF.vcf \
    -mode SNP \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file recal_GVCF.vcf \
    --tranches-file recal_GVCF.tranches \
    -O applied_recal_GVCF.vcf
echo "ApplyRecalibration PROCESSED"

#Annotate variants using snpEff"
java -jar $SNPEFF \
    -c $FILES_ANNOT/snpEff.config \
    -no-downstream \
    -no-upstream \
    -onlyProtein \
    Pf3D7v3 \
    applied_recal_GVCF.vcf > annotated.vcf   

#Include core regions from Pf genetic crosses version 1
bcftools annotate \
    -a $FILES_ANNOT/regions-20130225.onebased.txt.gz \
    -h $FILES_ANNOT/regions.hdr \
    -Ov \
    -o coreregions.vcf \
    -c CHROM,FROM,TO,RegionType annotated.vcf

#Annotate global barcode SNPs from Neafsey et al., 2008"
bcftools annotate \
    -a $FILES_ANNOT/global_barcode_tidy.txt.gz \
    -h $FILES_ANNOT/global_barcode.hdr \
    -Ov \
    -o barcodecoreregions.vcf \
    -c CHROM,FROM,TO,GlobalBarcode coreregions.vcf

#Select only biallelic SNPs
gatk SelectVariants \
    -R $REF_GEN \
    -V barcodecoreregions.vcf \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC \
    -O SNPs.vcf

echo "SelectVariants PROCESSED"

#Annotate VCF file with additional filters at the variant level
gatk VariantFiltration \
    -R $REF_GEN \
    --filter-name LowQualVQ -filter "VQSLOD <= 0.0" \
    --filter-name NotCore -filter "RegionType != 'Core'" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -V SNPs.vcf \
    -O SNP_filtered.vcf


#Clean VCF for analysis ready
gatk SelectVariants \
    -R $REF_GEN \
    -V SNP_filtered.vcf \
    -O SNP_filtered2.vcf \
    -select 'vc.isNotFiltered()'

#Exclude all genotypes with a filter flag not equal to "."(a missing value) or PASS
#vcftools --vcf $PATH_VCF/SNP_filtered2.vcf --remove-filtered-all --recode --stdout > $PATH_VCF/pass_SNPs.vcf
#echo "FIltering PROCESSED"

#SNP with a minor allele frequency of 3%
#vcftools --vcf SNP_filtered2.vcf --recode --recode-INFO-all --maf 0.03 --out SNP_filtered3.vcf

#missingness of at most 10% across the isolates
#bcftools view -i 'F_MISSING<0.1' $PATH_VCF/SNP_filtered3.vcf.recode.vcf > $PATH_VCF/SNP_filtered4.vcf
#bcftools view -i 'MIN(FMT/DP)>10' $PATH_VCF/SNP_filtered4.vcf.recode.vcf > $PATH_VCF/finalVCF.vcf


