#!/bin/bash

# This script should be run from a folder containing the a file named assemblies_suffix.txt containing 3 columns
#(1: name used in delta filenames (e.g. jagger), 2: abbreviated name for first line of chromosome level assemblies (e.g. jag), 3: full name (e.g. Jagger.1.1.))

#-----
# read in chromosome level delta files and unzip them (nucmer)
mkdir nucmerfiles#wget -np -r -P ./nucmerfiles/ -A '*chinese*.delta.gz' -R '*cadenza*,*claire*,*paragon*,*robigus*,*weebil*' https://opendata.earlham.ac.uk/wheat/under_license/toronto/Brinton_etal_2020-05-20-Haplotypes-for-wheat-breeding/nucmer/

# IF ONLY SELECTING CHINESE REFERENCE, change to '*chinese_v*', or IF ONLY FOR QUERY, '*v_chinese*'

#-----
# read in genome level assemblies and unzip them (wgas)
mkdir wgas
wget -np -r -P ./wgas/ -A 'Triticum_aestivum*pseudomolecules.fasta.gz' https://webblast.ipk-gatersleben.de/downloads/wheat/pseudomolecules/

#-----
# split these genome level assemblies by chromosome (clas)
declare -a CHROMNAMES=("chr1A" "chr1B" "chr1D" "chr2A" "chr2B" "chr2D" "chr3A" "chr3B" "chr3D" "chr4A" "chr4B" "chr4D" "chr5A" "chr5B" "chr5D" "chr6A" "chr6B" "chr6D" "chr7A" "chr7B" "chr7D" "chrUn")
for GZ in "./wgas/webblast.ipk-gatersleben.de/downloads/wheat/pseudomolecules/Triticum_aestivum"*"pseudomolecules"*;
do
	gunzip $GZ
done

mkdir clas
for REF_NAME in "./wgas/webblast.ipk-gatersleben.de/downloads/wheat/pseudomolecules/Triticum_aestivum"*"pseudomolecules"*;
do
	wga=${REF_NAME##*/}
	for REF_CHROM in "${CHROMNAMES[@]}";
	do
		sleep .01
		samtools faidx $REF_NAME $REF_CHROM > './clas/'$wga'.'$REF_CHROM'.fa'
	done
done

#-----
# show-snps

# First need to work out which 2 lines we're working with, and which chromosome
# Then need to replace the first line of the delta file to point to the two corresponding chromosome level assemblies (nucmer needs this)
# Then run show snps

mkdir snps
OLDIFS=$IFS

# Locate assemblies_suffix file (see top of script)
ASSEMBLIES_SUFFIX="assemblies_suffix.txt"

mkdir TestFolder

declare -a CHROMNUMS=("1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D")

for chrom in "${CHROMNUMS[@]}";
do
	chromfolder=./nucmerfiles/opendata.earlham.ac.uk/wheat/under_license/toronto/Brinton_etal_2020-05-20-Haplotypes-for-wheat-breeding/nucmer/"${chrom}"/
	counter=0

	for ZP in "${chromfolder}"*;
	do
		gunzip "${ZP}"
	done

	for files in "${chromfolder}"*;
	do
		counter=$((++counter))

		# Working out which lines are the reference and query sequences
		file=${files##*/}
		IFS='_.'
		if [[ "$file" == *"mattis"* ]];
		then
			read -ra splitfilename <<<"$file"
			if [[ "${splitfilename[0]}" == "sy" ]];
			then
				firstassembly="sy_mattis"
				secondassembly=${splitfilename[3]}
			else
				firstassembly=${splitfilename[0]}
				secondassembly="sy_mattis"
			fi
		else
			read -ra splitfilename <<<"$file"
			firstassembly=${splitfilename[0]}
			secondassembly=${splitfilename[2]}
		fi

		# Need to restore the IFS value.
		IFS=$OLDIFS

		# need to change first line of chromosome assemblies to the form >chr1A__ari
		FIRSTSUFFIX=$(grep $firstassembly $ASSEMBLIES_SUFFIX | cut -f 2)
		FIRST_REF=">chr"$chrom"__"$FIRSTSUFFIX

		SECONDSUFFIX=$(grep $secondassembly $ASSEMBLIES_SUFFIX | cut -f 2)
		SECOND_REF=">chr"$chrom"__"$SECONDSUFFIX

		# full name of lines used in assembly filenames
		FIRSTLONG=$(grep $firstassembly $ASSEMBLIES_SUFFIX | cut -f 3)
		SECONDLONG=$(grep $secondassembly $ASSEMBLIES_SUFFIX | cut -f 3)

		#paths to chromosome assemblies as those containing first and second
		firstpath=./clas/Triticum_aestivum_"${FIRSTLONG}"pseudomolecules.fasta.chr"${chrom}".fa
		secondpath=./clas/Triticum_aestivum_"${SECONDLONG}"pseudomolecules.fasta.chr"${chrom}".fa

		# use awk to rewrite the first line of nucmer files to be firstpath secondpath
		awk -v firstp="$firstpath" -v secondp="$secondpath" '{if(NR==1) {print firstp " " secondp} else {print $0}}' $files > ./TestFolder/"${file}".delta && mv ./TestFolder/"${file}".delta $files

		# use awk to rewrite the first line of clas files to be FIRST_REF or SECOND_REF
		awk -v firstref="$FIRST_REF" '{if(NR==1) {print firstref} else {print $0}}' $firstpath > ./TestFolder/${FIRST_REF}.fa && mv ./TestFolder/${FIRST_REF}.fa $firstpath
		awk -v secondref="$SECOND_REF" '{if(NR==1) {print secondref} else {print $0}}' $secondpath > ./TestFolder/${SECOND_REF}.fa && mv ./TestFolder/${SECOND_REF}.fa $secondpath

		# run show-snps
		if [ $(( $counter % 3 )) -ne 0 ];
		then
			show-snps -Clr -H -T $files > ./snps/${file}.snps &
		else
			show-snps -Clr -H -T $files > ./snps/${file}.snps
		fi

	done
done
