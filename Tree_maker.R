library(devtools)
library(tinytex)
library(dplyr)
library(phytools)
library(tidyverse)
library(compbio4all)
library(ggmsa)
library(rentrez)
library(seqinr)
library(ape)
library(msa)
library(Biostrings)
library(plotrix)

mycols <- colors()[c(8, 12, 16, 20, 24, 28, 32, 36:500)]

fasta_cleaner <- function(fasta_object, parse = TRUE){
  
  fasta_object <- sub("^(>)(.*?)(\\n)(.*)(\\n\\n)","\\4",fasta_object)
  fasta_object <- gsub("\n", "", fasta_object)
  
  if(parse == TRUE){
    fasta_object <- stringr::str_split(fasta_object,
                                       pattern = "",
                                       simplify = FALSE)
  }
  return(fasta_object[[1]])
}
entrez_fetch_list <- function(db, id, rettype, ...){
  
  #setup list for storing output
  n.seq <- length(id)
  list.output <- as.list(rep(NA, n.seq))
  names(list.output) <- id
  
  # get output
  for(i in 1:length(id)){
    list.output[[i]] <- rentrez::entrez_fetch(db = db,
                                              id = id[i],
                                              rettype = rettype)
  }
  
  
  return(list.output)
}
big_seq <- read.csv("Data/ur_table.csv", header = T) %>% dplyr::rename("16S" = X16S,
                                                                       "18S" = X18S,
                                                                       "ITS" = ITS.rDNA) 
big_seq <- select(big_seq, "Sample", "Canals", "Zooxanthellate", "Substrate", "Encrustment","Skeleton", "COI", "Colonial", "Macrocnemic.Brachycnemic" ,"16S", "18S", "ITS")
BS_COI <- na.omit(select(big_seq, "Sample", "COI"))
#BS_COI <- BS_COI %>% filter(Substrate != "Unknown")
BS_COI <- select(BS_COI, "Sample", "COI") 
BS_COI_names <- BS_COI$"Sample"
BS_COI$"ID" <- 1:nrow(BS_COI)

BS_COI_ID <- BS_COI$ID
BS_COI <- BS_COI$"COI"
BS_COI_list <- data.frame(c(1))

for(i in 1:length(BS_COI))
  BS_COI_list[[i]] <- fasta_cleaner(BS_COI[[i]], parse = F)
BSCOIl_vector <- rep(NA, length(BS_COI_list))
length(BSCOIl_vector)
for(i in 1:length(BSCOIl_vector)){
  BSCOIl_vector[i] <- BS_COI_list[[i]]
}

BSCOIl_vector_ID <- Biostrings::AAStringSet(BSCOIl_vector)


names(BSCOIl_vector_ID) <- BS_COI_ID



BSCOIl_align_ID <- msa(BSCOIl_vector_ID,
                       method = "Muscle")

BSCOIl_align_seqinr_ID <- msaConvert(BSCOIl_align_ID, 
                                     type = "seqinr::alignment")

BSCOIl_dist <- seqinr::dist.alignment(BSCOIl_align_seqinr_ID, 
                                      matrix = "identity")
BSCOIl_tree <- nj(BSCOIl_dist)


keyvector <- as.integer(BSCOIl_tree$tip.label)


#This unshuffles the names
BS_COI_names_shuffled <- c(1:length(BS_COI_names))
for(i in 1:length(BS_COI_names)){
  BS_COI_names_shuffled[i] <- BS_COI_names[keyvector[i]]
  
}

BSCOIl_tree$tip.label <- BS_COI_names_shuffled #Renames ID labels to specimen names


par(mar = c(0, 0, 0, 0))
plot.phylo(BSCOIl_tree, 
           type = "phylogram", 
           label.offset = 0.4, 
           cex = 0.6,
           font = 1,
           use.edge.length = F,
           )
#edgelabels(cex = 0.8)j

