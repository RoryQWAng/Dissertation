library(devtools)
library(tinytex)
library(dplyr)
library(phytools)
library(tidyverse)
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
big_seq <- read.csv("Data/big_seq2.csv", header = T) %>% dplyr::rename("16S" = X16S,
                                                                       "18S" = X18S,
                                                                       "ITS" = ITS.rDNA) 
big_seq <- select(big_seq, "Sample", "Canals", "Zooxanthellate", "Substrate", "Encrustment","Skeleton", "COI", "Colonial", "Macrocnemic.Brachycnemic" ,"16S", "18S", "ITS", "Ocean")
BS_COI <- na.omit(select(big_seq, "Sample", "COI", "Ocean"))
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
                       method = "ClustalW")

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



plot.phylo(BSCOIl_tree, 
           type = "phylogram", 
           use.edge.length = F)


Ocean_mode_og <- read.csv("Data/zo_character.csv", row.names=NULL) %>% dplyr::rename("16S" = X16S,
                                                                                     "18S" = X18S,
                                                                                     "ITS" = ITS.rDNA) 

Ocean_mode <- Ocean_mode_og[!is.na(Ocean_mode_og$"COI"),]
#crust_mode
#Ocean_mode$"Ocean"[is.na(Ocean_mode$"Ocean")] <- "Unknown"
Ocean_mode <- Ocean_mode[!is.na(Ocean_mode$"Ocean"),] 

#good_tip <- which(crust_mode$"Encrustment" != "Unknown")
Ocean_mode <- Ocean_mode %>% filter(Ocean != "Unknown")
#crust_mode
Ocean_mode <- setNames(Ocean_mode[,4],rownames(Ocean_mode))

Ocean_mode_shuffled <- c(1:length(Ocean_mode))
for(i in 1:length(BS_COI_names)){
  Ocean_mode_shuffled[i] <- Ocean_mode[keyvector[i]]
  
}
cbind(Ocean_mode_shuffled, BS_COI_names_shuffled)

Ocean_mode <- as.factor(Ocean_mode_shuffled)

plotTree(BSCOIl_tree,type="fan",fsize=0.8,ftype="i",lwd=2)
#length(levels(crust_mode))
#levels(crust_mode)
cols <- setNames(mycols,levels(Ocean_mode))[c(1:length(levels(Ocean_mode)))]
tiplabels(pie=to.matrix(Ocean_mode,
                        levels(Ocean_mode)),piecol=cols,cex=0.2)
to.matrix(Ocean_mode,
          levels(Ocean_mode))
add.simmap.legend(colors=cols,prompt=FALSE,x=1.5*par()$usr[3],
                  y=0*par()$usr[3],fsize=0.5)

fixed_tree <- root(BSCOIl_tree, outgroup = "Antipathozoanthus macaronesicus")
fixed_tree <- multi2di(fixed_tree)
is.rooted(fixed_tree)
is.ultrametric(fixed_tree)
is.binary.tree(fixed_tree)


branch_lengths <- branching.times(fixed_tree)
print(branch_lengths)
absolute_tree <- compute.brlen(fixed_tree, mode = "absolute") #Getting rid of all the negative branch lengths by changing them to their absolute value

fitER<-ace(Ocean_mode,absolute_tree,model="ER",type="discrete")
fitER
plot.phylo(absolute_tree, type="phylogram", label.offset = 0.03, cex = 0.5)
tiplabels(pie=to.matrix(Ocean_mode, levels(Ocean_mode)),piecol=cols,cex=0.1)
nodelabels(node=1:absolute_tree$Nnode+Ntip(absolute_tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.3, label.offset = 3)

add.simmap.legend(colors=cols,prompt=FALSE,x=1*par()$usr[1],
                  y=40*par()$usr[2.5],fsize=0.4)

