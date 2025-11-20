# Read the skyline output files

source("Rscripts/BaseScripts.R")
require(data.table)

# Compare the added areas of each peptide vs. 'Protein Abundance' outputs from skyline
dat1<-fread("Data/Overview.csv")
dat2<-fread("Data/MSStats_Input_Ofab_2025.csv")

colnames(dat2)[1]<-"Protein"
pro1.1<-dat1[Protein == "XP_020600399.1"]
pro1.2<-dat2[Protein == "XP_020600399.1"]


pro1.1[`Protein Result`== "XP_020600399.1 20250311_UWPRLumos_Richmond_Ofav1"]
# duplicated
#         Protein             Quantification Protein Abundance
#1: XP_020600399.1 Normalized Area: 6.4799E+5         9.5223E+6
#2: XP_020600399.1 Normalized Area: 4.1132E+6         9.5223E+6
#3: XP_020600399.1 Normalized Area: 6.4799E+5         9.5223E+6
#4: XP_020600399.1 Normalized Area: 4.1132E+6         9.5223E+6


total<-aggregate(pro1.2$Area, by=list(pro1.2$`File Name`), sum) 
total[total$Group.1=="20250311_UWPRLumos_Richmond_Ofav1.mzML",]
# 20250311_UWPRLumos_Richmond_Ofav1.mzML 39595664
# This includes duplicates


# Normalized Area = add all areas of fragment ions except precurosrs for each peptide
# Protein Abundance = Add normalized areas of all peptides for the protein   


#################################
# remove duplicates

dat1d<-dat1[!duplicated(dat1),] #1570616 rows -> 1073079 rows
dat2d<-dat2[!duplicated(dat2),] #13,561,915 -> 13139431
13561915-13139431
#422,484

pro1.2<-pro1.2[!duplicated(pro1.2),] #184 rows
pro1.1<-pro1.1[!duplicated(pro1.1),] #76 to 34?

#1. Ofav1 sample protein quantity from Overview
pro1.1[`Protein Result`== "XP_020600399.1 20250311_UWPRLumos_Richmond_Ofav1",c(1,5,10)]
#          Protein             Quantification Protein Abundance
#1: XP_020600399.1 Normalized Area: 6.4799E+5         9.5223E+6
#2: XP_020600399.1 Normalized Area: 4.1132E+6         9.5223E+6


#2. Ofav1 protein quantity from Emma's output
pro1.2[pro1.2$`File Name`=="20250311_UWPRLumos_Richmond_Ofav1.mzML",c(1,2,4,11,12)]
#          Protein Peptide Modified Sequence    Fragment Ion Normalized Area     Area
#1: XP_020600399.1                 VIPSLTSLK       precursor          647990        0
#2: XP_020600399.1                 VIPSLTSLK precursor [M+1]          647990  1374234
#3: XP_020600399.1                 VIPSLTSLK precursor [M+2]          647990        0
#4: XP_020600399.1                 VIPSLTSLK              y7          647990   331617
#5: XP_020600399.1                 VIPSLTSLK              y6          647990    94023
#6: XP_020600399.1                 VIPSLTSLK              y5          647990        0
#7: XP_020600399.1                 VIPSLTSLK              y4          647990        0
#8: XP_020600399.1                 VIPSLTSLK              y7          647990   222347
#9: XP_020600399.1                 AEGVVGIAR       precursor         4113200 10004489
#10: XP_020600399.1                 AEGVVGIAR precursor [M+1]         4113200  2834132
#11: XP_020600399.1                 AEGVVGIAR precursor [M+2]         4113200   823832
#12: XP_020600399.1                 AEGVVGIAR              y7         4113200  1316759
#13: XP_020600399.1                 AEGVVGIAR              y6         4113200        0
#14: XP_020600399.1                 AEGVVGIAR              y5         4113200  1080124
#15: XP_020600399.1                 AEGVVGIAR              y4         4113200   793741
#16: XP_020600399.1                 AEGVVGIAR              b3         4113200   922534

pro1.2_Ofav1<-pro1.2[pro1.2$`File Name`=="20250311_UWPRLumos_Richmond_Ofav1.mzML",]
# 2 peptides  VIPSLTSLK & AEGVVGIAR 
# normalized quantity 647990 4113200
sum(unique(pro1.2_Ofav1$`Normalized Area`))
# 4761190


## Add all normalized area to obtain protein abundance
pro1.2u<-pro1.2[!duplicated(pro1.2[,c("File Name","Peptide Modified Sequence")]),]

#replace NA with 0
pro1.2u$`Normalized Area`[is.na(pro1.2u$`Normalized Area`)]<-0


total1<-aggregate(pro1.2u$`Normalized Area`, by=list(pro1.2u$`File Name`), sum) #19 samples

total1


# How many peptides have NaN in Q value?

#check if the q value is the same for a given peptdie
pepdup<-dat2[!duplicated(dat2[,c("Peptide Modified Sequence","Detection Q Value")]),] #54267 unique peptides
pepdup2<-pepdup[duplicated(pepdup[,"Peptide Modified Sequence"]),]
# no duplicated peptides

nrow(pepdup[pepdup$`Detection Q Value`=="NaN",])
#24786

# Keep peptides with NaN for now.  



## We can use the ʻNormalized Areaʻ value for each peptide -> Each protein should have the sum of it. 
#replace NA with 
pep2$`Normalized Area`[is.na(pep2$`Normalized Area`)]<-0

proteinq<-aggregate(pep2$`Normalized Area`, by=list(pep2$`File Name`, pep2$Protein), sum)


# In Data-Independent Acquisition (DIA) proteomics analysis, DIA-NN software assigns statistical significance to identified precursors using q-values. A q-value of "NaN" (Not a Number) indicates that DIA-NN was unable to calculate a meaningful q-value for that particular precursor within the specified parameters

#In summary, a "NaN" q-value in DIA-NN indicates the software couldn't assign a reliable statistical significance to that precursor, likely due to data quality issues, spectral library limitations, or inappropriate data processing parameters. By examining the data and optimizing the analysis workflow, you can often address these issues and improve the overall identification and quantification performance. 

#######
# Overview2.csv -> following the advise, we removed the duplicated peptides

data<-fread("Data/Overview2.csv")
colnames(data)

#protein1<-data[Protein=="XP_020600399.1"]
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

colnames(data2)[colnames(data2)=="Protein Abundance Transition Averaged"]<-'Abundance'
data2$Abundance[is.na(data2$Abundance)]<-0
write.csv(data2, "Data/Ofav_protein_abundance.csv",row.names = F)

# calculate median for technical reps 
#Ofav1, Ofav23, Ofav43, Ofav60 have technical replicates

corals<-c("Ofav1", "Ofav26", "Ofav43", "Ofav60")
proteins<-data.frame(Protein=unique(data2$Protein))
Tech_median<-data.frame()
for (i in 1:length(corals)){
    df<-data2[data2$Sample==corals[i],]
    runs<-unique(df$Run)
    dfm<-df[df$Run==runs[1],]
    
    for(j in 1:nrow(dfm)){
        dfm$Abundance[j]<-median(df$Abundance[df$Protein==dfm$Protein[j]])
    }
    Tech_median<-rbind(Tech_median,dfm)
}
 
# Replace the tech rep individual data with median data

data3<-data2[!(data2$Sample %in% corals),]
data3<-rbind(data3, Tech_median)

write.csv(data3, "Data/Ofav_protein_abundance_techReps_medianed.csv",row.names = F)


# Run Paired comparison between two locations
# CFR vs. TR
locations<-c("CFR","MK","TR")


# Create a data frame with mean protein abundance for each location 

for (i in 1:nrow(proteins)){
    p<-proteins$Protein[i]
    df<-data3[data3$Protein==p,]
    proteins$CFR[i]<-mean(df$Abundance[df$Location==locations[1]])
    proteins$MK[i]<-mean(df$Abundance[df$Location==locations[2]])
    proteins$TR[i]<-mean(df$Abundance[df$Location==locations[3]])
    proteins$CFR_sd[i]<-sd(df$Abundance[df$Location==locations[1]])
    proteins$MK_sd[i]<-sd(df$Abundance[df$Location==locations[2]])
    proteins$TR_sd[i]<-sd(df$Abundance[df$Location==locations[3]])
}    

# calculate the fold change
for (i in 1:nrow(proteins[i])){
    proteins$CFR.MK.log2fold[i]<-log2(proteins$CFR[i]/proteins$MK[i])
    proteins$CFR.TR.log2fold[i]<-log2(proteins$CFR[i]/proteins$TR[i])
    proteins$MK.TR.log2fold[i]<-log2(proteins$MK[i]/proteins$TR[i])
    
}

changes<-proteins[,c("Protein","CFR.MK.log2fold","CFR.TR.log2fold","MK.TR.log2fold")]

# Manually calculated protein areas for each location
write.csv(proteins, "Output/Ofav.protein.areas.location.csv",row.names = F)


#CFK vs. TR

DF<-changes[,c(1,3)]
ggplor(DF, aes)
DF$abs.fold<-abs[]

for (i in 1:nrow(proteins)){
    p<-proteins$Protein[i]
    df<-data3[data3$Protein==p,]
    re1<-wilcox.test(df$Abundance[df$Location==locations[1]], df$Abundance[df$Location==locations[2]])
    re2<-wilcox.test(df$Abundance[df$Location==locations[1]], df$Abundance[df$Location==locations[3]])
    re3<-wilcox.test(df$Abundance[df$Location==locations[2]], df$Abundance[df$Location==locations[3]])
    
    t.test(df$Abundance[df$Location==locations[2]], df$Abundance[df$Location==locations[3]])
    
    df<-df[order(df$Location),]
    df$Sample<-factor(df$Sample, levels=df$Sample)
    ggplot(df, aes(x=Sample, y=Abundance, color=Location))+
        geom_point()
        
}

