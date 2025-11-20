library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)

get_symbol_from_name <- function(gene_name) {
    # Use HGNC search on the Name field (partial match)
    url <- paste0("https://rest.genenames.org/search/name/", URLencode(gene_name))
    r <- GET(url, add_headers(Accept = "application/json"))
    if (http_error(r)) return(tibble(name = gene_name, symbol = NA_character_))
    
    data <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = TRUE)
    
    # Turn docs into a tibble (if empty, return NA)
    docs <- tibble::as_tibble(data$response$docs)
    if (nrow(docs) == 0) return(tibble(name = gene_name, symbol = NA_character_))
    
    # Prefer Approved status and highest score (if _score exists)
    if (!"_score" %in% names(docs)) docs$`_score` <- NA_real_
    if (!"status" %in% names(docs)) docs$status <- NA_character_
    docs <- docs %>% arrange(desc(status == "Approved"), desc(`_score`))
    
    tibble(name = gene_name, symbol = docs$symbol[1])
}

map_names_to_symbols <- function(names_vec) {
    # be polite to the API and fail gracefully
    safe_get <- possibly(get_symbol_from_name, otherwise = tibble(name = NA, symbol = NA))
    map_dfr(names_vec, ~{
        Sys.sleep(0.05)
        safe_get(.x)
    })
}


# Example
gene_names <- c("tumor protein p53", "breast cancer type 1 susceptibility protein", "epidermal growth factor receptor")
result <- map_names_to_symbols(gene_names)
print(result)

ac<-read.csv("Data/Acerv.blast.table_updated.csv")

gene_names<-ac$Description
genes<-gsub("\\[.*?\\]",'', gene_names)
# remove -like
genes<-gsub("\\-like",'', genes)
# remove 'Predicted
genes<-gsub("PREDICTED\\:",'', genes)

genes<-trimws(genes)


result <- map_names_to_symbols(genes)
print(result)


# 




#Look at 'hypothetical protein"
hp<-result[grepl("hypothetical protein",result$name),]
# somehow it is filled with GPRASP1 
# Replace hypothetical protein with NA
result2<-result
result2$symbol[grepl("hypothetical protein",result2$name)]<-NA

# Add a gene symbol column
ac$GeneSymbol<-result2$symbol


add<-read.csv("Output/acerv.annotation_added.csv")
ac$GeneSymbol<-add$GeneSymbol
write.csv(ac, "Output/acerv.annotation.csv", row.names = F)



#UniprotIDs
uniprot<-read.table("Data/acerv.swissprot.uniprotID.txt")
colnames(uniprot)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

df2<-uniprot[!duplicated(uniprot$qseqid),] 

df2$UniprotID<-gsub("\\.\\d+",'',df2$sseqid)
df3<-df2[,c(1,3,13)]
colnames(df3)[1]<-"ProteinID"
ac<-ac[,c(1:7)]
ac2<-merge(ac, df3, by="ProteinID",all.x=TRUE)

write.csv(ac2,"Data/Acerv.blast.uniprot.csv", row.names = F)

hypo<-ac2[grepl("hypoth",ac2$Description),]




## Ofav blast results
uni<-fread("Data/ofav_swissprot_annot.txt")
colnames(uni)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")


prot<-read.csv("Data/Ofav_protein_list.csv")
prot[prot$Protein=="XP_020600396.1",]
uni2<-uni[!duplicated(uni$qseqid),] 

uni2$UniprotID<-gsub("\\.\\d+",'',uni2$sseqid)



uni3<-uni2[,c(1,3,13)]
colnames(uni3)[1]<-"ProteinID"
ofav_prot<-merge(prot, uni3, by.x="Protein", by.y="ProteinID",all.x=TRUE)
ofav_prot<-ofav_prot[!grepl("sp\\|",ofav_prot$Protein),]
range(ofav_prot$pident, na.rm = T)

write.csv(ofav_prot,"Output/Ofav.blast.uniprot.csv", row.names = F)

