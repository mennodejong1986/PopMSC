#!/bin/bash
# Script to detect haploblocks and to calculate raw genetic distances for each haploblock

# Depends on:
# - VCF_calcdist.perhaploblock.sh
# - VCF_phasing_splitdata.sh
# - VCF_calcdist.sh

# Usage:
# dos2unix VCF_haploblockanalysis.sh
# chmod +x VCF_haploblockanalysis.sh
# ./VCF_haploblockanalysis.sh


#####################################
BCFTOOLS=/opt/software/bcftools-1.16/bcftools
VCFTOOLS=/opt/software/vcftools/vcftools_0.1.17/bin/vcftools
SHAPEIT=/opt/software/shapeit4/shapeit4-4.2.2/bin/shapeit4.2
PLINK=/home/mdejong/software/plinkv2020-09-21/plink

MYCHROMS=mychroms3.txt		# head -37 /home/mdejong/bearproject/refgenome/brownbear_chrom/ASM358476v1_HiC.fasta.fai | cut -f1 > mychroms.txt
MYVCF=/home/mdejong/bearproject/snpfiles_february2022/Brown135.allsites.globalfilter.vcf.bgz
MYVCFSNPS=/home/mdejong/bearproject/snpfiles_february2022/Brown135.mysnps.missfilter.vcf.gz

indexvcf=FALSE
findblocks=FALSE
calcdistance=TRUE
#####################################



if [[ "$indexvcf" = TRUE ]]
	then
	echo "Indexing vcf files..."
	$BCFTOOLS index $MYVCF &
	$BCFTOOLS index $MYVCFSNPS $
	wait
	echo "Finished indexing."
fi


if [[ "$findblocks" = TRUE ]]
	then
	echo "Detecting haploblocks per chromosome..."
	for mychrom in $(cat $MYCHROMS)
		do
		echo $mychrom
		$BCFTOOLS view -r $mychrom $MYVCFSNPS -q 0.2:minor -O z -o mychrom.${mychrom}.vcf.gz
		$PLINK --vcf mychrom.${mychrom}.vcf.gz --allow-extra-chr --blocks no-pheno-req --out mychrom.${mychrom}.plink --blocks-max-kb 1500
		awk '$4>=25 && $5>50' mychrom.${mychrom}.plink.blocks.det | sed 's/ \+ / /g' | sed 's/ /\t/g' | cut -f1-5 | tail -n +2 > mychrom.${mychrom}.haploblocks.txt
		done
	wait
	echo "Haploblocks stored in files ending on 'haploblocks.txt'."
fi


if [[ "$calcdistance" = TRUE ]]
	then
	MYDIR=$(pwd)
	echo $MYDIR
	echo "Calculating pairwise raw genetic distances..."
	for mychrom in $(cat $MYCHROMS)
                do
		echo $mychrom
		if [ -d "$mychrom" ]
			then
			echo "Directory already exists."
			else	
			mkdir $mychrom
		fi
		cp mychrom.${mychrom}.haploblocks.txt $mychrom
		cp VCF_calcdist.sh $mychrom
		cp VCF_calcdist.perhaploblock.sh $mychrom
		cp VCF_phasing_splitdata.sh $mychrom
		#
		cd $mychrom
		MYSUBDIR=$(pwd)
		echo $MYSUBDIR
		./VCF_calcdist.perhaploblock.sh $mychrom & 
		cd $MYDIR
		done
fi


