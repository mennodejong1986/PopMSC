#!/bin/bash

#####################################
BCFTOOLS=/opt/software/bcftools-1.16/bcftools
MYVCF=/home/mdejong/bearproject/snpfiles_february2022/Brown135.mysnps.missfilter.vcf.gz

mychrom=$1
#####################################

echo $mychrom

nblocks=$(wc -l mychrom.${mychrom}.haploblocks.txt | cut -f1 -d ' ')
echo "chrom startpos endpos maf0 maf0.02 maf0.05 maf0.1 maf0.2" > mychrom.${mychrom}.haploblocks.nsites.txt

for blocknr in $(seq 1 $nblocks)
#for blocknr in $(seq 1 2)
	do
    	echo "Starting analyses on "$blocknr" out of "$nblocks" haploblocks in total."
    	blockstart=$(awk -v myline="$blocknr" 'NR==myline' mychrom.${mychrom}.haploblocks.txt | cut -f2)
    	blockend=$(awk -v myline="$blocknr" 'NR==myline' mychrom.${mychrom}.haploblocks.txt | cut -f3)
    	echo $blockstart" - "$blockend
	nsnvs0=$($BCFTOOLS view -r ${mychrom}:${blockstart}-${blockend} --exclude-types indels $MYVCF | grep -v '#' | wc -l)
	nsnvs1=$($BCFTOOLS view -r ${mychrom}:${blockstart}-${blockend} -q 0.02:minor --exclude-types indels $MYVCF | grep -v '#' | wc -l)
	nsnvs2=$($BCFTOOLS view -r ${mychrom}:${blockstart}-${blockend} -q 0.05:minor --exclude-types indels $MYVCF | grep -v '#' | wc -l)
	nsnvs3=$($BCFTOOLS view -r ${mychrom}:${blockstart}-${blockend} -q 0.1:minor --exclude-types indels $MYVCF | grep -v '#' | wc -l)
	nsnvs4=$($BCFTOOLS view -r ${mychrom}:${blockstart}-${blockend} -q 0.2:minor --exclude-types indels $MYVCF | grep -v '#' | wc -l)
	echo "$mychrom $blockstart $blockend $nsnvs0 $nsnvs1 $nsnvs2 $nsnvs3 $nsnvs4" >> mychrom.${mychrom}.haploblocks.nsites.txt
	done
echo "Done" $mychrom
