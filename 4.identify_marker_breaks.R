###use Rscript identify_marker_breaks.R ancestry-probs_rec_X_v5.txt ancestryinfer Dtei TUZ35_CT03

arrArgs<-commandArgs(trailingOnly = TRUE);
genos<-as.character(arrArgs[1])
path<-as.character(arrArgs[2])
species<-as.character(arrArgs[3])
cross<-as.character(arrArgs[4])

if(length(arrArgs)<2){
stop("usage is: Rscript identify_marker_breaks.R genotypes_file_name path_to:transpose_nameout.pl\n");
}#print usage

#transpose, not of R (but run thorugh it)
#perl scripts_mareyMap/transpose_nameout.pl ancestry-probs_rec_X_v5.txt

command1<-paste("perl ",path,"/transpose_nameout.pl ",genos,sep="")
system(command1)
#command1

#get the name of the transposed geno file
#ancestry-probs_rec_X_v5.txt_transposed
#genos_transposed<-"ancestry-probs_rec_X_v5.txt_transposed"

genos_transposed<-paste(genos,"_transposed",sep="")
#genos_transposed

#select markers from transposed file
#cut -f 1 ancestry-probs_rec_X_v5.txt_transposed | perl -p -e 's/:/\t/g' | cut -f 1 | uniq | tail -n +2 > ancestry-probs_rec_X_v5.txt_header

select_markers<-paste("cut -f 1 ",genos_transposed," | perl -p -e 's/:/\t/g' | cut -f 1 | uniq | tail -n +2 > ",genos,"_header",sep="")
system(select_markers)
#select_markers

#get chromosomes
#chroms_raw<-read.csv("ancestry-probs_rec_X_v5.txt_header",sep="\t",head=FALSE)
chroms_raw<-read.csv(file=paste(genos,"_header",sep=""),sep="\t",head=FALSE)

chroms<-chroms_raw$V1
#chroms

#create outfile name, if outfile exists already in the same dir it will be replaced
#outfile<-paste("ancestry-probs_rec_X_v5.txt_ancestrytransitions_allchrs")

outfile<-paste(genos,"_ancestrytransitions_allchrs",sep="")
file.remove(outfile)

for(k in 1:length(chroms)){

command2=paste("grep ",chroms[k]," ",genos_transposed," | perl -p -e ","'","s/:/\t/g","'"," > ",genos_transposed,"_",chroms[k],sep="")
system(command2)
#command2

file=paste(genos_transposed,"_",chroms[k],sep="")
#file

data<-read.csv(file=file,sep="\t",head=FALSE,as.is=T)
#head(data)

intervals<-{}
count=0

for(y in 3:(length(data[1,]))){

focal_geno=data[,y]
sites<-data[,2]
na_count=0;
geno_prev=subset(focal_geno,!is.na(focal_geno)==TRUE)[1];
first_region=1
last_geno=focal_geno[1]
start=sites[1]

for(x in 1:length(focal_geno)){
geno_current=focal_geno[x];

if(!is.na(geno_current)==TRUE){

	if(geno_current == geno_prev){

	} else{

	if(is.na(last_geno)==TRUE){
	stop=sites[x]
	geno_prev=geno_current

	intervals<-rbind(intervals,cbind(start,stop,y-2,geno_prev))

	}#make sure previous genotype was NA for this mode

	}#geno_current does not equal geno_previous

} else{

if((is.na(geno_current)==TRUE)	&(is.na(last_geno)==FALSE)){
	start=sites[x-1]
}#option 1: this is the first site in the interval region

}#if is/is not NA


if((is.na(geno_current)==FALSE)&(is.na(last_geno)==FALSE)&(geno_current != last_geno)){
	start=sites[x-1]
	stop=sites[x]
	geno_prev=geno_current

	intervals<-rbind(intervals,cbind(start,stop,y-2,geno_prev))
	count=count+1
}#option 2: the interval was between this marker and the last marker


last_geno=geno_current

}#for all sites in the focal individual

}#for all individuals

write.table(cbind(rep(as.character(chroms[k]),length(intervals[,1])),intervals),file=outfile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)

command3<-paste(genos,"_",chroms[k],sep="")

}#all chroms


#final format

to_format_DF <- read.table (outfile, header =F)

#because there are 0s this doesnt work 
#to_format_DF$previous_geno <- ifelse(to_format_DF$V5 == 1, 2, 1)

#format_DF <- to_format_DF[, c(1, 2, 3, 4, 6,5)]
format_DF <- to_format_DF[, c(1, 2, 3, 4,5)]

colnames(format_DF)<-c("Chromossome","marker1","marker2","sample_nr","marker2_geno")

#get the actual sample names, instead of just interval numbers
command4=paste("awk '{print $1}' ",genos," | sed 's/id/sample_name/g' > ",genos,"_list",sep="")
system(command4)

#read the list
individual_names_list<-paste(genos,"_list",sep="")
individual_names_list_R <- read.table (individual_names_list, header =T)

#add sample numbers
individual_names_list_R$sample_nr<-1:nrow(individual_names_list_R)
colnames (individual_names_list_R)<-c("sample_name","sample_nr")

#merge
final_table_geno_transitions <- merge (format_DF,individual_names_list_R, by="sample_nr", all=TRUE)
final_table_geno_transitions$Species <-species
final_table_geno_transitions$Cross <-cross

#final_table_geno_transitions
final_table_geno_transitions_format <- final_table_geno_transitions[, c(7,8,2, 3, 4, 5, 1, 6)]

write.table(as.data.frame(final_table_geno_transitions_format), file=outfile,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


#clean enviroment
#remove header file
file_header=paste(genos,"_header",sep="")
file.remove(file_header)

#remove genos transposed
file.remove(genos_transposed)

#remove individual names list
file.remove(individual_names_list)
