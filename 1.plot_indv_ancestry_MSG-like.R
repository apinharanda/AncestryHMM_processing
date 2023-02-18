####v2 run from command line
#AP Feb 2022
#this plots the ancestry in an MSG-like away 
#it plots it for each sample that we want
#also outputs a table with genomes by sample at each position. This can be used to compare parent calls to AIMs (to correct genotypes)

#conda R activate


##how to run it

#Rscript plot_indv_ancestry_MSG-like.R DteiTUZ35_ancestry-probs_0.5_all X _read_1_val_1.fq _read_2_val_4.fq par2_L5_TUZ35_10 GenotypeTrack 


rm(list=ls())

library(data.table)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)

arrArgs<-commandArgs(trailingOnly = TRUE);

file_tsv_rec <- as.character(arrArgs[1]) ###the file is ancestry-probs_2L_2R_3L_3R_X.tsv_rec.txt
file_tsv_rec

Chr <- as.character(arrArgs[2]) #add chr name
Chr

read1end <- as.character(arrArgs[3]) #add the end of read 1 (check name of individuals in rec file), something like _read_1_val_1.fq
read1end

read2end <- as.character(arrArgs[4]) #add the end of read 1 (check name of individuals in rec file), something like _read_2_val_4.fq
read2end

sampleID <- as.character(arrArgs[5]) #add sample ID
sampleID

out_name <- as.character(arrArgs[6]) #name of plot
out_name


##to test
#X_025_09 <- fread("DteiTUZ35_ancestry-probs_0.5_all",header=T) 
#Chr<-"X"
#read1end<-"_read_1_val_1.fq"
#read2end<-"_read_2_val_4.fq"
#sampleID<-"par2_L5_TUZ35_10"
#out_name<-"GenotypeTrack"


##import table
X_025_09 <- fread(file_tsv_rec,header=T)


# Transpose table
X_025_09.T <- t(X_025_09[,2:ncol(X_025_09)])

# Set the column headings from the first column in the original table
head_df<-as.data.frame(X_025_09[,1] )
t_X_025_09_df <- as.data.frame(X_025_09.T)
colnames(t_X_025_09_df) <- head_df [,1]

#change names
for_str =paste(Chr,":\\d{1}", sep="")

names(t_X_025_09_df) <- gsub(read1end, "", names(t_X_025_09_df), fixed = TRUE)
names(t_X_025_09_df) <- gsub(read2end, "", names(t_X_025_09_df), fixed = TRUE)

#select chr
t_X_025_09_df_selected<-t_X_025_09_df[str_detect(row.names(t_X_025_09_df), for_str), ]

for_rows = paste(Chr,":", sep="")
rownames(t_X_025_09_df_selected) <- gsub(for_rows, "", rownames(t_X_025_09_df_selected), fixed = TRUE)

t_X_025_09_df_selected$pos <- rownames(t_X_025_09_df_selected)
t_X_025_09_df_selected$pos_num<-as.numeric(as.character(t_X_025_09_df_selected$pos))
t_X_025_09_df_selected$pos <- NULL


t_X_025_09_df_long <-
t_X_025_09_df_selected %>%
  pivot_longer(!pos_num, names_to = "Sample", values_to = "Genotype")


# creating a column to dataframe based on values in other column:

t_X_025_09_df_long <- t_X_025_09_df_long %>% 
  mutate(plot_geno_par1_het = if_else(Genotype == 1, "0", "NA"))

t_X_025_09_df_long <- t_X_025_09_df_long %>% 
  mutate(plot_geno_par2_het = if_else(Genotype == 1, "2", "NA"))


t_X_025_09_df_long <- t_X_025_09_df_long %>% 
  mutate(plot_geno_par2_homo = if_else(Genotype == 2, "2", "NA"))

t_X_025_09_df_long <- t_X_025_09_df_long %>% 
  mutate(plot_geno_par1_homo = if_else(Genotype == 0, "0", "NA"))

##
t_X_025_09_df_long$Genotype<-NULL

t_X_025_09_df_long_2 <-
t_X_025_09_df_long %>%
  pivot_longer(cols = -c(pos_num,Sample), names_to = "what", values_to = "Genotype")

noNA_table<-t_X_025_09_df_long_2[!grepl("NA", t_X_025_09_df_long_2$Genotype),]
noNA_table$Genotype_num <-as.numeric(as.character(noNA_table$Genotype))



full_out=paste(out_name, Chr,sampleID,".pdf", sep="_")
title= paste(out_name, Chr, sampleID, sep="_")



pdf (full_out)
ggplot()+
geom_hline(aes(yintercept = 1),colour = "grey",size=2)+
geom_hline(aes(yintercept = 1),colour = "black")+

#these were to add annotations for the simulation checks
#annotate("text", y = 2.2, x = 12474642/1000000, label = "*",colour="mediumseagreen",fontface=2,size=7)+
#geom_vline(aes(xintercept =  12474642/1000000),colour="mediumseagreen",linetype="dotted")+
#annotate("text", y = 2.2, x = 1791325/1000000, label = "*",colour="mediumseagreen",fontface=2,size=7)+
#geom_vline(aes(xintercept = 1791325/1000000),colour="mediumseagreen",linetype="dotted")+
#annotate("text", y = 2.2, x = 5602912/1000000, label = "*",colour="mediumseagreen",fontface=2,size=7)+
#geom_vline(aes(xintercept = 5602912/1000000),colour="mediumseagreen",linetype="dotted")+

geom_rect(data=subset(noNA_table, Sample %in% c(sampleID)& what %in% c("plot_geno_par2_homo")),size=0.3,aes(xmin=(pos_num/1000000),xmax=(pos_num/1000000),ymin=Genotype_num-1, ymax=Genotype_num,colour=as.character(Genotype),fill=as.character(Genotype)))+
geom_rect(data=subset(noNA_table, Sample %in% c(sampleID)& what %in% c("plot_geno_par1_homo")),size=0.3,aes(xmin=(pos_num/1000000),xmax=(pos_num/1000000),ymin=Genotype_num+1,ymax=Genotype_num,colour=as.character(Genotype),fill=as.character(Genotype)))+
geom_point(data=subset(noNA_table, Sample %in% c(sampleID)& what %in% c("plot_geno_par2_het")),shape = 108,size=0.3,aes(x=(pos_num/1000000),y=Genotype_num,colour=as.character(Genotype)))+
geom_point(data=subset(noNA_table, Sample %in% c(sampleID)& what %in% c("plot_geno_par1_het")),shape = 108,size=0.3,aes(x=(pos_num/1000000),y=Genotype_num,colour=as.character(Genotype)))+
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
labs(x ="Position (Mb)", y = "Genotype")+
scale_fill_manual(name = "Genotype", values = c("0" = "red", "2" = "blue"), labels = c("Par1", "Par2"))+
scale_colour_manual(name = "", values = c("0" = "red", "2" = "blue"),labels = c("", "", "",""))+
scale_y_continuous(name="Genotypes", limits=c(0,2.25),breaks=c(0,1,2),labels=c("par1", "","par2"))+
ggtitle(title)
dev.off()


##output table for AIMs cross ref
#this may be helpful when masking the genomes

#name_out = paste(out_name, Chr,sampleID,".genotypes_at_eachPos.tsv", sep="_")

#write.table(noNA_table, file=name_out, sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)


