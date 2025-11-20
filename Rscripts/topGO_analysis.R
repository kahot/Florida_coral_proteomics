# GO

library(topGO)

# inputs you prepare:
# 1) all_genes: named numeric/binary vector over the *background* (1=hit, 0=not)
#    names(all_genes) are Acerv gene IDs present in your experiment
# 2) gene2GO: list mapping gene_id -> character vector of GO terms

acerv_prots<-read.csv("Output/acerv.annotation_added.csv")

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = all_genes,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO,
              nodeSize = 10)