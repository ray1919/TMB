#!/bin/bash
# Author: zzqr@163.com
# Date: 2020-03-12
# calculate TMB inspired by TMB.py

# dependency
# - bedtools
# - bcftools
jannovar=/home/zhaorui/miniconda3/bin/jannovar
jannovar_data_hg19=/home/zhaorui/bg904/tools/jannovar/data/hg19_refseq.ser
TSG=/home/zhaorui/bg904/tools/tumor_analysis/TMB/src/TSG.txt
somatic=/home/zhaorui/bg904/tools/tumor_analysis/TMB/src/somatic.tsv
germline_common=/home/zhaorui/bg904/tools/tumor_analysis/TMB/src/germline_common.tsv
exome_hg19_bed=/home/zhaorui/bg904/tools/tumor_analysis/TMB/src/Twist_Exome_Target_hg19.bed.gz

# temporary file produced during processing
# - biallelic.vcf
# - biallelic.ex.vcf
# - biallelic.jv.vcf
# - TSG-stop_gained.tsv
# - mutation.tsv
# - mutation_50.tsv

# stop on error
set -e

# Usage: ./TMB_calculate.sh SmallVariants.vcf target_region.bed

SmallVariants=$1
target_region=$2

## expand collapsed vcf
bcftools norm -m - $SmallVariants -o biallelic.vcf 2>/dev/null

# intersect with exome region
bedtools intersect -header -a biallelic.vcf -b $exome_hg19_bed > biallelic.ex.vcf

## annotate
$jannovar annotate-vcf --report-no-progress -d $jannovar_data_hg19 -i biallelic.ex.vcf -o biallelic.jv.vcf 2>/dev/null

paste <(awk '$7=="PASS"' biallelic.jv.vcf | cut -f1,2,4,5) \
      <(awk '$7=="PASS"' biallelic.jv.vcf | cut -f8 | awk 'BEGIN{FS="|"}{printf "%s\t%s\n",$2,$4}') | \
      grep stop_gained | \
      grep -w -f $TSG | \
      cut -f1-4 | \
      sed 's/^chr//' > TSG-stop_gained.tsv

## Mutation
paste <(awk '$7=="PASS"' biallelic.ex.vcf | cut -f1,2,4,5) \
      <(awk '$7=="PASS"' biallelic.ex.vcf | cut -f8 | awk 'BEGIN{FS=";"}{print $1}' | awk 'BEGIN{FS="="}{print $2}') | \
      awk '$5>=0.05' | cut -f1-4 | sed 's/^chr//' > mutation.tsv

# Mutation 50
paste <(awk '$7=="PASS"' biallelic.ex.vcf | cut -f1,2,4,5) \
      <(awk '$7=="PASS"' biallelic.ex.vcf | cut -f8 | awk 'BEGIN{FS=";"}{print $1}' | awk 'BEGIN{FS="="}{print $2}') | \
      awk '$5<0.5 && $5>=0.05' | cut -f1-4 | sed 's/^chr//' > mutation_50.tsv

# count of valid mutations in TMB
cTMB=`grep -wf mutation.tsv $germline_common | \
    grep -v -wf - mutation.tsv | \
    grep -v -wf $somatic | \
    grep -v -wf TSG-stop_gained.tsv | \
    wc -l`

# count of valid mutations (AF < 0.5) in TMB
cTMB_50=`grep -wf mutation_50.tsv $germline_common | \
    grep -v -wf - mutation_50.tsv | \
    grep -v -wf $somatic | \
    grep -v -wf TSG-stop_gained.tsv | \
    wc -l`

# panel size in Mb
pSIZE=`bedtools merge -i $target_region | bedtools intersect -a - -b $exome_hg19_bed | awk '{len=$3-$2;sum+=len}END{print sum/1000000}'`

# TMB_50
TMB_50=`echo $cTMB_50/$pSIZE | bc -l`

# TMB
TMB=`echo $cTMB/$pSIZE | bc -l`
echo $TMB $TMB_50 | awk '{if($2 > 20){print $1}else{print $2}}'
