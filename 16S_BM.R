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
BS_16S <- na.omit(select(big_seq, "Sample", "16S", "Macrocnemic.Brachycnemic"))
BS_16S <- BS_16S %>% filter(Macrocnemic.Brachycnemic != "Unknown")
BS_16S <- BS_16S %>% filter(Macrocnemic.Brachycnemic != "Unkown")
BS_16S <- select(BS_16S, "Sample", "16S") 
BS_16S_names <- BS_16S$"Sample"
BS_16S$"ID" <- 1:nrow(BS_16S)

BS_16S_ID <- BS_16S$ID
BS_16S <- BS_16S$"16S"
BS_16S_list <- data.frame(c(1))

for(i in 1:length(BS_16S))
  BS_16S_list[[i]] <- fasta_cleaner(BS_16S[[i]], parse = F)
BS16Sl_vector <- rep(NA, length(BS_16S_list))
length(BS16Sl_vector)
for(i in 1:length(BS16Sl_vector)){
  BS16Sl_vector[i] <- BS_16S_list[[i]]
}

BS16Sl_vector_ID <- Biostrings::AAStringSet(BS16Sl_vector)


names(BS16Sl_vector_ID) <- BS_16S_ID



BS16Sl_align_ID <- msa(BS16Sl_vector_ID,
                       method = "ClustalW")

BS16Sl_align_seqinr_ID <- msaConvert(BS16Sl_align_ID, 
                                     type = "seqinr::alignment")

BS16Sl_dist <- seqinr::dist.alignment(BS16Sl_align_seqinr_ID, 
                                      matrix = "identity")
BS16Sl_tree <- njs(BS16Sl_dist)


keyvector <- as.integer(BS16Sl_tree$tip.label)


#This unshuffles the names
BS_16S_names_shuffled <- c(1:length(BS_16S_names))
for(i in 1:length(BS_16S_names)){
  BS_16S_names_shuffled[i] <- BS_16S_names[keyvector[i]]
  
}

BS16Sl_tree$tip.label <- BS_16S_names_shuffled #Renames ID labels to specimen names



plot.phylo(BS16Sl_tree, 
           type = "phylogram", 
           use.edge.length = F)


BM_og <- read.csv("Data/zo_character.csv", row.names=NULL) %>% dplyr::rename("16S" = X16S,
                                                                                   "18S" = X18S,
                                                                                   "ITS" = ITS.rDNA) 

BM_mode <- BM_og[!is.na(BM_og$"16S"),]
#crust_mode
#BM_mode$"BM"[is.na(BM_mode$"BM")] <- "Unknown"
BM_mode <- BM_mode[!is.na(BM_mode$"Macrocnemic.Brachycnemic"),] 

#good_tip <- which(crust_mode$"Encrustment" != "Unknown")
BM_mode <- BM_mode %>% filter(Macrocnemic.Brachycnemic != "Unknown")
BM_mode <- BM_mode %>% filter(Macrocnemic.Brachycnemic != "Unkown")
#crust_mode
BM_mode <- setNames(BM_mode[,13],rownames(BM_mode))

BM_mode_shuffled <- c(1:length(BM_mode))
for(i in 1:length(BS_16S_names)){
  BM_mode_shuffled[i] <- BM_mode[keyvector[i]]
  
}
cbind(BM_mode_shuffled, BS_16S_names_shuffled)

BM_mode <- as.factor(BM_mode_shuffled)

plotTree(BS16Sl_tree,type="fan",fsize=0.8,ftype="i",lwd=2)
#length(levels(crust_mode))
#levels(crust_mode)
cols <- setNames(mycols,levels(BM_mode))[c(1:length(levels(BM_mode)))]
tiplabels(pie=to.matrix(BM_mode,
                        levels(BM_mode)),piecol=cols,cex=0.2)
to.matrix(BM_mode,
          levels(BM_mode))
add.simmap.legend(colors=cols,prompt=FALSE,x=1.5*par()$usr[3],
                  y=0*par()$usr[3],fsize=0.5)

fixed_tree <- root(BS16Sl_tree, outgroup = "Antipathozoanthus macaronesicus")
fixed_tree <- multi2di(fixed_tree)
is.rooted(fixed_tree)
is.ultrametric(fixed_tree)
is.binary.tree(fixed_tree)


branch_lengths <- branching.times(fixed_tree)
print(branch_lengths)
absolute_tree <- compute.brlen(fixed_tree, mode = "absolute") #Getting rid of all the negative branch lengths by changing them to their absolute value

fitER<-ace(BM_mode,absolute_tree,model="ER",type="discrete")
fitER
plot.phylo(absolute_tree, type="phylogram", label.offset = 0.03, cex = 0.5)
tiplabels(pie=to.matrix(BM_mode, levels(BM_mode)),piecol=cols,cex=0.3)
nodelabels(node=1:absolute_tree$Nnode+Ntip(absolute_tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.3, label.offset = 3)

add.simmap.legend(colors=cols,prompt=FALSE,x=1*par()$usr[1],
                  y=3.2*par()$usr[2.5],fsize=0.8)

