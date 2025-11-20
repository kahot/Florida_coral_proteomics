# Read the skyline output files

source("Rscripts/BaseScripts.R")
require(data.table)
library(ggplot2)

#dat2<-fread("Data/MSStats_Input_Ofab_2025.csv")

#######
# Overview2.csv -> following the advise, we removed the duplicated peptides in skyline

data<-fread("Data/Overview2.csv")
colnames(data)

# protein description table
prot<-data[!duplicated(data$Protein),]
prot<-prot[,c("Protein","Protein Description")]
write.csv(prot,"Data/Ofav_protein_list.csv", row.names = F)

protein1<-data[Protein=="XP_020600399.1"]
#each row corresponds to each peptide. 
# "Protein Abundance Transition Averaged" is the same as "protein Abundance"

# Format the data for MSStats analysis:
data$`Normalized Area`[is.na(data$`Normalized Area`)]<-0

meta<-read.csv('Data/20250311_UWPRLumos_Richmond_sequence.csv')
data<-merge(data, meta[1:2], by.x = 'Replicate',by.y = 'File.Name', all.x = TRUE)

# Remove space and a period from Sample IDs:
data$Sample<-gsub('\\.','',data$Sample.ID)
data$Sample<-gsub(' ','',data$Sample)

# Add "annotation" information required by MSstats (Condition, BioReplicate, Run)
sampleinfo<-read.csv('Data/Florida_proteins.csv')
sampleinfo$Sample<-paste0('Ofav',sampleinfo$ID)
colnames(sampleinfo)[colnames(sampleinfo)=="Collection.Location"]<-'Location'


#Ofav54 must of Ofav59? 
data$Sample[data$Sample=="Ofav54"]<-"Ofav59"
data<-merge(data,sampleinfo[,c('Location','Sample')],by='Sample', all.x=TRUE )

# Create run number (unique to each File.Name)
data$Run<-as.numeric(factor(data$`File Name`))

# Make one row per protein and keep only the protein abundance info.
#number of unqiue proteins 
length(unique(data$Protein)) #9669

# number of samples
length(unique(data$Sample)) #15

# number of runs
length(unique(data$Run)) #19
# Some samples have technical duplicates


data2<-data[,c('Sample','Protein','Peptide','Protein Description','Protein Abundance Transition Averaged','Location','File Name','Run')]


data2<-data2 %>%
    distinct(Run, Protein, .keep_all = TRUE)
# nrow(data_filtered)/19
#183711 rows

colnames(data2)[colnames(data2)=="Protein Abundance Transition Averaged"]<-'Abundance'


#data2<-data2[order(data2$Location),]
#data2$Run<-factor(data2$Run, levels=unique(data2$Run))
#ggplot(data2, aes(x=Run, y=Abundance, fill=Location))+
#    geom_bar(stat="identity",width=0.8)

# Check which one has technical replicates
samples<-data2[,c("Sample",'Run')]
samples<-samples[!duplicated(samples$Run),]

#Lok at Ofav1 (Run 1 & 8)
ofav1<-data2[data2$Sample=="Ofav1",]
ofav1.df<-data.frame(Protein=unique(ofav1$Protein))
run1<-ofav1[ofav1$Run==1]
run8<-ofav1[ofav1$Run==8]
ofav1.df<-merge(ofav1.df, run1[,c("Protein","Abundance")], by="Protein")
colnames(ofav1.df)[colnames(ofav1.df)=="Abundance"]<-"Run1"

ofav1.df<-merge(ofav1.df, run8[,c("Protein","Abundance")], by="Protein")
colnames(ofav1.df)[colnames(ofav1.df)=="Abundance"]<-"Run8"
ofav1.df$difference<-ofav1.df$Run1-ofav1.df$Run8
ggplot(ofav1.df[1:500,],aes(x=Protein, y=difference ))+
    geom_point()

ofav1m<-melt(ofav1.df[,1:3])
remove(raw2)
