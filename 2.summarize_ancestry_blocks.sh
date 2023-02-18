#!/bin/bash

# v1 AP Jan 4 2021
#script to output size of UNFILTERED ancestry blocks from Ancestry-HMM output after running Rscript scripts_mareyMap/identify_marker_breaks.R 
#need to know how many columns there are in file 

#doesnt necessary need to be run, its just a summary of the data

#so, for example do this first
#awk -F'\t' '{print NF; exit}'  ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt_transposed_X
#1154


# e.g.
# bash summarize_ancestry_blocks.sh 1 1154 ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt_transposed_X ancestry_blocks_X_Dtei_unfiltered.txt


##### read in command line

#first column to summarise (where the first sample geno appears)
sample_1=$1 #I cant seem to be able to add it as numeric, just leave this part as hard coded
#until when in the file to summarise samples (which is the column with the last sample to analyse)
sample_last=$2  #I cant seem to be able to add it as numeric, just leave this part as hard coded
#file input
in_file=$3
#output file where the ancestry blocks UNFiltered will be written to
out_file=$4



#separate file by columns file transposed after running first script, ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt_transposed_X
rm -rf *tmp*

########ADD number of columns here

for i in {3..194}; do
    awk -v x=$i '{print $1"\t"$2"\t" $x }' $in_file | grep -v "NA" > sample_${i}.tmp.all
done

#get the pos where the sample changes genos
#the first and last are missing
##also add number of columns here
for i in {3..194}; do
    awk -v x=$i '{print $1"\t"$2"\t" $x }' $in_file | grep -v "NA" | awk 'BEGIN  { getline; cmp=$3; preline=$0 } cmp!=$3{ print preline, $0 } { cmp=$3; preline=$0 }' OFS='\n' > sample_${i}.tmp.middle
done

#change file names to start at 1

i=1
for file in $(ls *.tmp.middle | sort -nk 1.8) ; do 
     mv -f $file $(echo $file | sed 's/[0-9]*.tmp.middle$/'${i}'.tmp.middle/')
     i=$((++i))
done

i=1
for file in $(ls *.tmp.all | sort -nk 1.8) ; do 
     mv -f $file $(echo $file | sed 's/[0-9]*.tmp.all$/'${i}'.tmp.all/')
     i=$((++i))
done


#so get the first row
for file in $(ls *.tmp.all | sort -nk 1.8) ; do
    head -n1 "$file" >> all.tmp.top
done

#and get the last row so this can be appended to
for file in $(ls *.tmp.all | sort -nk 1.8) ; do
    tail -n1 "$file" >> all.tmp.bottom
done

#separate first and last row files for each sample by sample name, creating new files
awk '!/^$/{print  > "sample_"NR".all.tmp.top" }' all.tmp.top
awk '!/^$/{print  > "sample_"NR".all.tmp.bottom" }' all.tmp.bottom


#list the three sets of files
ls -1 sample_*middle | sort -nk 1.8 > middle_files_tmp
ls -1 sample_*top | sort -nk 1.8 > top_files_tmp
ls -1 sample_*bottom | sort -nk 1.8 > bottom_files_tmp

#for each sample cat the start, middle and end files
paste top_files_tmp middle_files_tmp bottom_files_tmp | awk '{print "cat " $1 " "$2" "$3" > "$1".middle.bottom"}' | sed '1 i#!/bin/bash' > make_indiv.files.tmp.sh

bash make_indiv.files.tmp.sh

#add sample name to each row
for f in *.all.tmp.top.middle.bottom
do
    fname=$(basename "$f")
    fname=${fname##filename_}
    echo "$fname"
    sed -i "s/\$/\t$fname/" "$fname"
done

#cat the several individual files
cat sample_*.all.tmp.top.middle.bottom | sed 's/.all.tmp.top.middle.bottom//g' > out.file.tmp

awk '{printf "%s%s",$0,(NR%2?FS:RS)}' out.file.tmp | awk '{print $1"\t"$2"\t"$6"\t"$3"\t"$7"\t"$4}' | sed '1 iChr\tBlock_start\tBlock_end\tGeno_blockStart\tGeno_blockEnd\tSample' > $out_file

#remove temporary files
rm -rf *tmp*


