###Feb 2022 04
#v2 AP run from the command line
#v 3 changed plotting Feb17

#this script calculates and plots the probability of homozygosity for par 2 scaffolds for each chr
#run first
#perl ./ancestryinfer/parsetsv_to_genotypes_v2.pl ancestry-probs-par1_2L_2R_3L_3R_X.tsv ancestry-probs-par2_2L_2R_3L_3R_X.tsv 0.9 ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt
#need to input the read ends
#and the chr accordingly
#
#Rscript prob_homo_par2.R DteiTUZ35_ancestry-probs_0.5_all X _read_1_val_1.fq _read_2_val_4.fq DteiTUZ_ancestryProbs

#conda R activate

rm(list=ls())

library(data.table)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)

arrArgs<-commandArgs(trailingOnly = TRUE);

file_tsv_rec <- as.character(arrArgs[1]) #this is the input file which is the output of the script ./ancestryinfer/parsetsv_to_genotypes_v2.pl 
file_tsv_rec
Chr <- as.character(arrArgs[2]) #add chr name
Chr
read1end <- as.character(arrArgs[3]) #add the end of read 1 (check name of individuals in rec file), something like _read_1_val_1.fq
read1end
read2end <- as.character(arrArgs[4]) #add the end of read 1 (check name of individuals in rec file), something like _read_2_val_4.fq
read2end
out_name <- as.character(arrArgs[5]) #name of plot
out_name

#direct input
#X_025_09 <- fread("ancestry-probs_2L_2R_3L_3R_X.tsv_rec_males.txt",header=T) 
#Chr <- "X"
#read1end <- "_read_1_val_1.fq"
#read2end <- "_read_2_val_4.fq"
#out_name <- "DteiTUZ_ancestryProbs_unfilteredmales"
 
 
#read table 
X_025_09 <- fread(file_tsv_rec,header=T) 

# Transpose table
X_025_09.T <- t(X_025_09[,2:ncol(X_025_09)])

# Set the column headings from the first column in the original table
head_df<-as.data.frame(X_025_09[,1] )
t_X_025_09_df <- as.data.frame(X_025_09.T)
colnames(t_X_025_09_df) <- head_df [,1]

#change names
names(t_X_025_09_df) <- gsub(read1end, "", names(t_X_025_09_df), fixed = TRUE)
names(t_X_025_09_df) <- gsub(read2end, "", names(t_X_025_09_df), fixed = TRUE)

#select chr
for_str =paste(Chr,":\\d{1}", sep="")

t_X_025_09_df_selected<-t_X_025_09_df[str_detect(row.names(t_X_025_09_df),for_str ), ]

for_rows = paste(Chr,":", sep="")
rownames(t_X_025_09_df_selected) <- gsub(for_rows, "", rownames(t_X_025_09_df_selected), fixed = TRUE)

t_X_025_09_df_selected$pos <- rownames(t_X_025_09_df_selected)
t_X_025_09_df_selected$pos_num<-as.numeric(as.character(t_X_025_09_df_selected$pos))
t_X_025_09_df_selected$pos <- NULL

t_X_025_09_df_long <-
t_X_025_09_df_selected %>%
  pivot_longer(!pos_num, names_to = "Sample", values_to = "Genotype")


t_X_025_09_df_long_count_genos <- t_X_025_09_df_long %>% group_by(pos_num) %>%
    count(Genotype)  
    
t_X_025_09_df_long_count_genos_noNA <- t_X_025_09_df_long_count_genos[!is.na(t_X_025_09_df_long_count_genos$Genotype), ]

t_X_025_09_df_long_count_genos_noNA_sum <-  t_X_025_09_df_long_count_genos_noNA %>% group_by(pos_num)%>% 
summarise(samples_with_genotypes = sum(n))

freq<- merge (t_X_025_09_df_long_count_genos_noNA, t_X_025_09_df_long_count_genos_noNA_sum, by ="pos_num")

freq$freq<-freq$n/freq$samples_with_genotypes


full_out=paste(out_name, Chr,".pdf", sep="_")
title= paste(out_name, Chr, "Prob par2 alleles", sep="\n")

pdf(full_out)
ggplot(data=freq) +
geom_hline(aes(yintercept =  0.5),colour="grey",linetype="dotted")+
#geom_point(aes(x=pos_num/1000000,y= freq, group=Genotype),size=.001,alpha=0.3)+ 
geom_line(aes(x=pos_num/1000000,y= freq, group=Genotype,colour=as.character(as.numeric(Genotype))),size=.2)+ 
#facet_wrap(~Genotype, ncol=1, scales="free")+
theme_bw()+
theme(aspect.ratio=0.3)+
theme(axis.line = element_line(colour = "black"),
aspect.ratio=0.3,
panel.border=element_rect(colour="black", size=0.5),
text=element_text(size=11),
legend.text = element_text(size = 10),
legend.title = element_text(face = "bold", size = 10),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
axis.title.y= element_text(vjust=1,face="bold",size = 10, margin=margin(10,20,10,0)),
axis.title.x = element_text(vjust=1,margin=margin(10,20,10,0),face="bold", size = 10),
axis.text.x = element_text(size=9),
axis.text.y = element_text(size=9),
strip.background = element_blank(),
legend.position="bottom",
strip.text.x = element_text(face = "bold"),
strip.text.y = element_text(face = "bold"),
plot.title = element_text(face="bold"))+
labs(x ="Chr position (Mb)", y = "Mean probability")+
scale_colour_manual(name="Genotype",values=c("red","purple","blue"))+
ggtitle(title)+
#xlim(0,24)+
ylim(0,1)+
dev.off() 

mean <- freq %>%
  group_by(Genotype) %>%
  summarise(meanfreq = mean(freq), sd = sd(freq))  
  
full_out=paste(out_name, Chr,".txt", sep="_")

write.table(as.data.frame(mean), file=full_out,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


