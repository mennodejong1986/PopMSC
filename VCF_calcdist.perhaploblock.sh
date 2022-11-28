#!/bin/bash
# Script to calculate raw genetic distances per haploblock

# Expects to find in working directory:
# - VCF_phasing_splitdata.bash
# - VCF_calcdist.sh

# Usage:
# ./VCF_calcdist.perhaploblock.sh inputdata.vcf.gz


#####################################
BCFTOOLS=/opt/software/bcftools-1.16/bcftools
VCFTOOLS=/opt/software/vcftools/vcftools_0.1.17/bin/vcftools
SHAPEIT=/opt/software/shapeit4/shapeit4-4.2.2/bin/shapeit4.2	# this script uses by default Beagle instead of Shapeit (because shapeit excludes multi-allelic during phasing and imputation)
BEAGLE=/home/mdejong/software/beagle/beagle.22Jul22.46e.jar	
PLINK=/home/mdejong/software/plinkv2020-09-21/plink

MYVCF=/home/mdejong/bearproject/snpfiles_february2022/Brown135.allsites.globalfilter.vcf.bgz
MYVCFSNPS=/home/mdejong/bearproject/snpfiles_february2022/Brown135.mysnps.filter.vcf.gz
MYREF=/home/mdejong/bearproject/refgenome/brownbear_chrom/ASM358476v1_HiC.fasta

use_shapeit=FALSE
ncores=4		# number of threads used by shapeit or beagle
nindsvcf=270

mychrom=$1
#####################################


# $BCFTOOLS index $MYVCFSNPS
# echo "Finding haploblocks..."
# $BCFTOOLS view -r ${mychrom} $MYVCFSNPS -q 0.2:minor -O z -o mychrom.${mychrom}.vcf.gz
# $PLINK --vcf mychrom.${mychrom}.vcf.gz --allow-extra-chr --blocks no-pheno-req --out mychrom.${mychrom}.plink --blocks-max-kb 1500
# awk '$4>=25 && $5>50' mychrom.${mychrom}.plink.blocks.det | sed 's/ \+ / /g' | sed 's/ /\t/g' | cut -f1-5 | tail -n +2 > mychrom.${mychrom}.haploblocks.txt

 
nblocks=$(wc -l mychrom.${mychrom}.haploblocks.txt | cut -f1 -d ' ')
if [ -f vcfdist.all.txt ]; then rm vcfdist.all.txt; fi
touch vcfdist.all.txt


for blocknr in $(seq 1 $nblocks)
	do
	echo "Starting analyses on "$blocknr" out of "$nblocks" haploblocks in total."
	blockstart=$(awk -v myline="$blocknr" 'NR==myline' mychrom.${mychrom}.haploblocks.txt | cut -f2)
	blockend=$(awk -v myline="$blocknr" 'NR==myline' mychrom.${mychrom}.haploblocks.txt | cut -f3)
	echo $blockstart" - "$blockend
	#
	# EXTRACT REGION:
	echo "Extracting region..."
	$BCFTOOLS view -r ${mychrom}:${blockstart}-${blockend} --exclude-types indels $MYVCF -O z > myhaploblock.vcf.gz
	echo "Indexing..."
	$BCFTOOLS index myhaploblock.vcf.gz
	#
	# SPLIT MULTI-ALLELIC SITES:
	# $BCFTOOLS norm --multiallelics -snps myhaploblock.vcf.gz -O z > myhaploblock.split.vcf.gz &
	# $BCFTOOLS index myhaploblock.split.vcf.gz
	#
	# PHASE DATA:
	echo "Phasing..."
	# whatshap phase -o mysnps.phased.whatshap.vcf.gz --reference=${MYREF} --chromosome ${mychrom} myhaploblock.split.vcf.gz in1.bam in2.bam in3.bam etc
	if [[ "$use_shapeit" = TRUE ]]
		then
		echo "...using shapeit (warning: biallelic snps only)"
		$SHAPEIT --input myhaploblock.vcf.gz --region ${mychrom} --output myhaploblock.phased.vcf.gz -T $ncores
		else
		echo "...using beagle"
		java -Xmx64g -jar $BEAGLE gt=myhaploblock.vcf.gz out=myhaploblock.phased.beagle nthreads=${ncores} chrom=${mychrom}
		# Note: unlike Shapeit4, Beagle output all sites (monomorphic and polymorphic)
		# To speed up distance calculations, we remove monomorphic sites.
		# This can later easily be corrected for, because missing data has been imputed, so all individuals have the same number of total sites.
		$BCFTOOLS index myhaploblock.phased.beagle.vcf.gz
		$BCFTOOLS view --threads 2 --min-alleles 2 -O z myhaploblock.phased.beagle.vcf.gz -o myhaploblock.phased.vcf.gz	
	fi
	# CONVERT FORMAT:
	# set flag 'haploid' to TRUE:
	echo "Converting data into haploid format..."
	./VCF_phasing_splitdata.sh
	#
	# CALCULATE RAW GENETIC DISTANCE
	# set flag 'haploiddata' to TRUE.
	echo "Calculating raw genetic distances..."
	./VCF_calcdist.sh 1 $nindsvcf
	#
	echo "Adding results to output file..."
	grep -v 'startbp' vcfdist.1_${nindsvcf}.txt > vcfdist.1_${nindsvcf}.noheader.txt
	cat vcfdist.all.txt vcfdist.1_${nindsvcf}.noheader.txt > vcfdist.all.tmp.txt
	mv vcfdist.all.tmp.txt vcfdist.all.txt
	echo "Finished region."
	done

cp vcfdist.all.txt vcfdist.${mychrom}.txt
echo "DONE. All analyses finished. All raw genetic distances stored in file 'vcfdist.mychrom.txt'."
