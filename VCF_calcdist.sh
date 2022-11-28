#!/bin/bash
# A shell script to calculate uncorrected pairwise distance between any two samples in a vcf file.

# USAGE:
# Execute by typing, for example:
# ./VCF_calcdist 1 10
# This will calculate all pairwise scores for individuals 1 to 10.
# For example, if the dataset contains in total 100 individuals, then the command above will calculate the scores for 10*100 = 1000 pairwise comparisons.

# If you forget to add the two numbers, you will receive the error:
# line 78: ((: i = : syntax error: operand expected (error token is "= ")


# Run the command in two steps. First convert the data, afterwards do the actual analysis.
# This is especially important when running multiple analyses at the same time (using 'VCF_calcdist_parallel.sh').

# Note that data will be converted to an unzipped geno file.
# This means that you should not use a big vcf file as input, otherwise the geno file will be gigantic.
# Instead, use for example a thinned vcf file.

# First use bcftools query to extract sample names and genotype scores from vcf file:
# /opt/software/bcftools/bcftools-1.9/bin/bcftools query --list-samples allsites.thinned.200.vcf.gz > myinput.samples.txt
# /opt/software/bcftools/bcftools-1.9/bin/bcftools query -f '[\t%GT]\n' allsites.thinned.200.vcf.gz | sed 's/^[ \t]*//' > myinput.geno.txt




##############################################################
start1=$1       	# Specify this number on the command line.
end1=$2         	# Specify this number on the command line.

start2=1
end2=270        	# Set this number to the total of individuals in the input vcf file.

convertdata=TRUE
run_loop=TRUE
haploiddata=TRUE
haplodiploiddata=FALSE

MYVCF=mydata.split.haploid.vcf.gz
BCFTOOLS=/opt/software/bcftools-1.16/bcftools
##############################################################




if [[ "$convertdata" == TRUE ]]
	then
	echo "Converting vcf to geno..."
	$BCFTOOLS query --list-samples $MYVCF > myinput.samples.txt
	$BCFTOOLS query -f '[\t%GT]\n' $MYVCF | sed 's/^[ \t]*//' > myinput.geno.txt
	zgrep -v '#' $MYVCF | head -1 | cut -f1-2 > myinput.pos.txt
	# example commands:
	# autosomal:
	# /opt/software/bcftools/bcftools-1.9/bin/bcftools query --list-samples allsites.thinned.200.vcf.gz > myinput.samples.txt
	# /opt/software/bcftools/bcftools-1.9/bin/bcftools query -f '[\t%GT]\n' allsites.thinned.200.vcf.gz | sed 's/^[ \t]*//' > myinput.geno.txt
	# Y-chromosome:
	# /opt/software/bcftools/bcftools-1.9/bin/bcftools query --list-samples allsites.globalfilter.vcf.gz > myinput.samples.txt &
	# /opt/software/bcftools/bcftools-1.9/bin/bcftools query -f '[\t%GT]\n' allsites.globalfilter.vcf.gz | sed 's/^[ \t]*//' > myinput.geno.txt &
	else
	echo "Assuming input files 'myinput.geno.txt' and 'myinput.samples.txt' are already present."
fi

if [[ "$run_loop" == TRUE ]]
	then
	if [[ "$haploiddata" == TRUE ]]
		then
		echo "ind1 ind2 name1 name2 nmiss n0 n1 chrom startbp" > vcfdist.${start1}_${end1}.txt
		else
		if [[ "$haplodiploiddata" == TRUE ]]
                	then
			echo "ind1 ind2 name1 name2 nmiss h0 h1 d0 d1 d2he d2ho hd0 hd1 hd2 chrom startbp" > vcfdist.${start1}_${end1}.txt
			else
			echo "ind1 ind2 name1 name2 nmiss n0 n1 n2he n2ho chrom startbp" > vcfdist.${start1}_${end1}.txt
		fi
	fi
	mychrom=$(cut -f1 myinput.pos.txt)
	mypos=$(cut -f2 myinput.pos.txt)
	echo "Starting pairwise comparisons..."
	for (( i = $start1; i <= $end1; i++ ))
	do
		echo "$i"
		ind1=$(awk -v myline="$i" 'NR==myline' myinput.samples.txt)
		cut -f${i} myinput.geno.txt > myind1.${i}.txt
		for (( j = $start2 ; j <= $end2; j++ ))
		do
			# select data:
			# echo "$i $j"
			ind2=$(awk -v myline="$j" 'NR==myline' myinput.samples.txt)
			cut -f${j} myinput.geno.txt > myind2.${i}_${j}.txt
			paste myind1.${i}.txt myind2.${i}_${j}.txt | sed 's|\/|\t|g' > mypair.${i}_${j}.txt
			grep -v '\.' mypair.${i}_${j}.txt > mypair.${i}_${j}.nomissing.txt
			nmiss=$(grep '\.' mypair.${i}_${j}.txt | wc -l)
			
			# count number of sites:
			if [[ "$haploiddata" == TRUE ]]
				then
				# haploid data:
				n0=$(awk '$1==$2' mypair.${i}_${j}.nomissing.txt | wc -l)
				n1=$(awk '$1!=$2' mypair.${i}_${j}.nomissing.txt | wc -l)
				# write output:
				echo "$i $j $ind1 $ind2 $nmiss $n0 $n1 $mychrom $mypos" >> vcfdist.${start1}_${end1}.txt
				else
				if [[ "$haplodiploiddata" == TRUE ]]
					then
					# haplodiploid data:
					head -1 myind1.${i}.txt > myind1.${i}.line1.txt
					head -1 myind2.${i}_${j}.txt > myind2.${i}_${j}.line1.txt
					indbool1=$(grep '/' myind1.${i}.line1.txt | wc -l)
					indbool2=$(grep '/' myind2.${i}_${j}.line1.txt | wc -l)
					if [[ "$indbool1" == 1 && "$indbool2" == 1  ]]
						then
						echo "Individual $i and $j are both diploid."
						h0=0
						h1=0
						d0=$(awk '$1==$3 && $2==$4 && $1==$2' mypair.${i}_${j}.nomissing.txt | wc -l)
                                        	d1=$(awk '($1!=$3 && $2==$4) || ($1==$3 && $2!=$4)' mypair.${i}_${j}.nomissing.txt | wc -l)
                                        	d2he=$(awk '$1==$3 && $2==$4 && $1!=$2' mypair.${i}_${j}.nomissing.txt | wc -l)
                                        	d2ho=$(awk '$1!=$3 && $2!=$4' mypair.${i}_${j}.nomissing.txt | wc -l)
						hd0=0
                                                hd1=0
                                                hd2=0
						echo "$i $j $ind1 $ind2 $nmiss $h0 $h1 $d0 $d1 $d2he $d2ho $hd0 $hd1 $hd2 $mychrom $mypos" >> vcfdist.${start1}_${end1}.txt
					fi		
					if [[ "$indbool1" == 1 && "$indbool2" == 0  ]]
						then
						echo "Individual $i is diploid and $j is haploid."
						h0=0
						h1=0
						d0=0
						d1=0
						d2he=0
						d2ho=0
						hd0=$(awk '$1==$3 && $2==$3' mypair.${i}_${j}.nomissing.txt | wc -l)
						hd1=$(awk '$1==$3 && $2!=$3 || $1!=$3 && $2==$3' mypair.${i}_${j}.nomissing.txt | wc -l)
						hd2=$(awk '$1!=$3 && $2!=$3' mypair.${i}_${j}.nomissing.txt | wc -l)
						echo "$i $j $ind1 $ind2 $nmiss $h0 $h1 $d0 $d1 $d2he $d2ho $hd0 $hd1 $hd2 $mychrom $mypos" >> vcfdist.${start1}_${end1}.txt				
					fi	
					if [[ "$indbool1" == 0 && "$indbool2" == 1  ]]
						then
						echo "Individual $i is haploid and $j is diploid."
						h0=0
                                                h1=0
                                                d0=0
                                                d1=0
                                                d2he=0
                                                d2ho=0
                                                hd0=$(awk '$1==$2 && $1==$3' mypair.${i}_${j}.nomissing.txt | wc -l)
                                                hd1=$(awk '$1==$2 && $1!=$3 || $1!=$2 && $1==$3' mypair.${i}_${j}.nomissing.txt | wc -l)
                                                hd2=$(awk '$1!=$2 && $1!=$3' mypair.${i}_${j}.nomissing.txt | wc -l)
                                                echo "$i $j $ind1 $ind2 $nmiss $h0 $h1 $d0 $d1 $d2he $d2ho $hd0 $hd1 $hd2 $mychrom $mypos" >> vcfdist.${start1}_${end1}.txt						
					fi	
					if [[ "$indbool1" == 0 && "$indbool2" == 0 ]]
                                                then
                                                echo "Individual $i and individual $j are both haploid."
                                        	h0=$(awk '$1==$2' mypair.${i}_${j}.nomissing.txt | wc -l)
                                		h1=$(awk '$1!=$2' mypair.${i}_${j}.nomissing.txt | wc -l)
						d0=0
						d1=0
						d2he=0
						d2ho=0
						hd0=0
						hd1=0
						hd2=0
						echo "$i $j $ind1 $ind2 $nmiss $h0 $h1 $d0 $d1 $d2he $d2ho $hd0 $hd1 $hd2 $mychrom $mypos" >> vcfdist.${start1}_${end1}.txt
					fi					
					rm myind1.${i}.line1.txt myind2.${i}_${j}.line1.txt
					else
					# diploid data:
					n0=$(awk '$1==$3 && $2==$4 && $1==$2' mypair.${i}_${j}.nomissing.txt | wc -l)
					n1=$(awk '($1!=$3 && $2==$4) || ($1==$3 && $2!=$4)' mypair.${i}_${j}.nomissing.txt | wc -l)
					n2he=$(awk '$1==$3 && $2==$4 && $1!=$2' mypair.${i}_${j}.nomissing.txt | wc -l)
					n2ho=$(awk '$1!=$3 && $2!=$4' mypair.${i}_${j}.nomissing.txt | wc -l)
					# write output:
					echo "$i $j $ind1 $ind2 $nmiss $n0 $n1 $n2he $n2ho $mychrom $mypos" >> vcfdist.${start1}_${end1}.txt
				fi	
			fi
			rm myind2.${i}_${j}.txt mypair.${i}_${j}.txt mypair.${i}_${j}.nomissing.txt
		done
		rm myind1.${i}.txt
	done
	echo "Finished analyses. Data stored in file 'vcfdist.txt'."
	sed -i 's/ /\t/g' vcfdist.${start1}_${end1}.txt
	echo "For diploid data, dxy can be calculated as follows: (n1*0.5+n2he*0.5+n2ho)/(n0+n1+n2he+n2ho)."
	echo "For haploid data, the formula is simply: n1/(n0+n1)."
	echo "For haplodiploid data, a combination of these formulas is required."
	echo "For any of these three data types, the calculation can be performed using the function 'add2inds2()' implemented in SambaR."
	echo "This Sambar function will automatically detect the data type and perform the calculations accordingly."
	echo "See also the script 'VCF_calcdist_plotinR.txt'."
	else
	echo "Not running analyses because the flag 'run_loop' is set to FALSE."
fi
