# AncestryHMM_processing
Scripts to plot the output of ancestryHMM
a series of scripts that allow to look at the output of AncestryHMM
these have been written specifically to analyse low coverage data (introgression project for example)

################################################################################################################################################################################################
################################################################################################################################################################################################

1.plot_indv_ancestry_MSG-like.R

#this plots the samples in an MSG like way
#run for the example data as like this

##how to run it
#
#Rscript plot_indv_ancestry_MSG-like.R DteiTUZ35_ancestry-probs_0.5_all X _read_1_val_1.fq _read_2_val_4.fq par2_L5_TUZ35_10 GenotypeTrack 

#the arguments are:
#1) the file is ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt
#2) Chr to plot
#3) the end of the sample name (usually read info) for read 1
#4) the end of the sample name (usually read info) for read 2, when there is just one read (ie one end) write the same thing twice 
#5) the sample ID to plot 
#6) any extra descriptor for the plot

#Specifically for the examples in this folder, after finishing ancestry HMM run
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_1 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_2 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_3 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_4 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_5 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_6 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_7 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_8 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_9 GenotypeTrack_simulation
Rscript plot_indv_ancestry_MSG-like.R ancestry-probs_allchrs.tsv_rec.txt X _AAAAAA _AAAAAA indivA6QDA15_10 GenotypeTrack_simulation

#the results are inside this folder

################################################################################################################################################################################################
################################################################################################################################################################################################

2.summarize_ancestry_blocks.sh

#Script to output size of UNFILTERED ancestry blocks from Ancestry-HMM
#This script uses the direct output from AncestryHMM and summarises ancestry blocks for each sample
#it outputs tables that can then be used to filter the data, it is not always necessary to run, its more of a way to check the data

#Because temp files are created while it is running only one chr can be run at a time (ie need to wait until one finishes to run another one)
#It is still a script under improvement - the number of columns of each file are hardcoded inside the script

#Need to know how many columns there are in file 
#so, for example do this first
awk -F'\t' '{print NF; exit}'  ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt_transposed_X
#1154

#then inside the script change the number of columns in two locations

#Improvements will include identifying the number of columns within the script (or adding it to the command line)

#to run use
bash summarize_ancestry_blocks.sh 3 1154 ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt_transposed_X ancestry_blocks_X_Dtei_unfiltered.txt

#where
##### read in command line
#1)first column to summarise (where the first sample geno appears)
#I cant seem to be able to add it as numeric, just leave this part as hard coded
#
#2)until when in the file to summarise samples (which is the column with the last sample to analyse)
#I cant seem to be able to add it as numeric, just leave this part as hard coded
#
#3)file in (from AncestryHMM output)
#
#4) output file where the ancestry blocks UNFiltered will be written to


################################################################################################################################################################################################
################################################################################################################################################################################################

3.prob_homo_par2.R

#this script plots ancestry probabilities at each position averaged by all individuals that could be genotyped
##need to

conda activate R

#the arguments are:
#1) the file is ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt
#2) Chr to plot
#3) the end of the sample name (usually read info) for read 1
#4) the end of the sample name (usually read info) for read 2, when there is just one read (ie one end) write the same thing twice 
#5) any extra descriptor for the plot

Rscript prob_homo_par2.R ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt 2L _read_1_val_1.fq _read_2_val_4.fq DteiTUZ_ancestryProbs &
Rscript prob_homo_par2.R ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt 2R _read_1_val_1.fq _read_2_val_4.fq DteiTUZ_ancestryProbs &
Rscript prob_homo_par2.R ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt 3L _read_1_val_1.fq _read_2_val_4.fq DteiTUZ_ancestryProbs &
Rscript prob_homo_par2.R ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt 3R _read_1_val_1.fq _read_2_val_4.fq DteiTUZ_ancestryProbs &
Rscript prob_homo_par2.R ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt X _read_1_val_1.fq _read_2_val_4.fq DteiTUZ_ancestryProbs &

#it is plotting both points and a line, sometimes both look good, sometimes one looks better than the other
#change script inside to test that if necessary


################################################################################################################################################################################################
################################################################################################################################################################################################

4.identify_marker_breaks.R

#identifies marker breaks to then be used to filter the recombination map
#needs to be run before 5

Rscript identify_marker_breaks_2.R ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt scripts_mareyMap Dsim SE

################################################################################################################################################################################################
################################################################################################################################################################################################

