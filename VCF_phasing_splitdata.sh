#!/bin/bash
# Script to split the haplotypes in a phased vcf-file into separate genotypes in haploid or diploid format 

# For haploid data: 0|1 becomes 0 1
# For diploid data: 0|1 becomes 0/0 1/1.
# Note: for diploid data, conversion works for biallelic data only.
# In case of multi-allelic data, convert with this script to haploid data, and subsequently use VCF_haploid2diploid script to convert to diploid.

# Samples names are appended with -1 and -2.    


#######################################
MYVCF=myhaploblock.phased.vcf.gz
haploid=TRUE				# should the data be converted to haploid (1, 0) or diploid data (1/1, 0/1, 0/0)?
#######################################

zgrep -v '#' $MYVCF > mydata.tmp.txt
zgrep '##' $MYVCF > myheader.tmp.txt
zgrep -v '##' $MYVCF| head -1 > myheaderline.tmp.txt

echo "Editing sample names..."
cut -f1-9 myheaderline.tmp.txt > myheaderline.first9columns.tmp.txt
cut -f10- myheaderline.tmp.txt | sed 's/\t/\n/g' | sed 's/$/-1/' > mysamples1.tmp.txt
cut -f10- myheaderline.tmp.txt | sed 's/\t/\n/g' | sed 's/$/-2/' > mysamples2.tmp.txt
paste -d '\n' mysamples1.tmp.txt mysamples2.tmp.txt | tr '\n' '\t' | sed 's/[ \t]*$//' > mysamples3.tmp.txt
paste myheaderline.first9columns.tmp.txt mysamples3.tmp.txt > myheaderline.new.tmp.txt 

echo "Editing data..."
if [[ "$haploid" = TRUE ]]
	then
	echo "Converting to haploid format..."
	# 24-09-2022: depreciated command suitable for biallelic data only:
	# sed 's/0|1/0>1/g' mydata.tmp.txt | sed 's/1|0/1>0/g' | sed 's/0|0/0>0/g' | sed 's/1|1/1>1/g' | sed 's/>/\t/g' > mydata.new.tmp.txt
	sed 's/|/\t/g' mydata.tmp.txt > mydata.new.tmp.txt
	echo "Combining..."
        cat myheader.tmp.txt myheaderline.new.tmp.txt mydata.new.tmp.txt | gzip > mydata.split.haploid.vcf.gz
	rm my*tmp.txt
	echo "Done :-) Output is stored in file called 'mydata.split.haploid.vcf.gz'."
	else
	echo "Converting to diploid format..."
	# 24-09-2022: command is suitable for biallelic data only:
	sed 's/0|1/0%0>1%1/g' mydata.tmp.txt | sed 's/1|0/1%1>0%0/g' | sed 's/0|0/0%0>0%0/g' | sed 's/1|1/1%1>1%1/g' | sed 's|%|/|g' | sed 's/>/\t/g' > mydata.new.tmp.txt
	echo "Combining..."
	cat myheader.tmp.txt myheaderline.new.tmp.txt mydata.new.tmp.txt | gzip > mydata.split.diploid.vcf.gz
	rm my*tmp.txt
	echo "Done :-) Output is stored in file called 'mydata.split.diploid.vcf.gz'."	
fi

