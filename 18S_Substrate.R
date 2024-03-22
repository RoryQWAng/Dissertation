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
BS_18S <- na.omit(select(big_seq, "Sample", "18S", "Substrate"))
BS_18S <- BS_18S %>% filter(Substrate != "Unknown")
BS_18S <- select(BS_18S, "Sample", "18S") 
BS_18S_names <- BS_18S$"Sample"
BS_18S$"ID" <- 1:nrow(BS_18S)

BS_18S_ID <- BS_18S$ID
BS_18S <- BS_18S$"18S"
BS_18S_list <- data.frame(c(1))

for(i in 1:length(BS_18S))
  BS_18S_list[[i]] <- fasta_cleaner(BS_18S[[i]], parse = F)
BS18Sl_vector <- rep(NA, length(BS_18S_list))
length(BS18Sl_vector)
for(i in 1:length(BS18Sl_vector)){
  BS18Sl_vector[i] <- BS_18S_list[[i]]
}

BS18Sl_vector_ID <- Biostrings::AAStringSet(BS18Sl_vector)


names(BS18Sl_vector_ID) <- BS_18S_ID



BS18Sl_align_ID <- msa(BS18Sl_vector_ID,
                       method = "ClustalW")

BS18Sl_align_seqinr_ID <- msaConvert(BS18Sl_align_ID, 
                                     type = "seqinr::alignment")

BS18Sl_dist <- seqinr::dist.alignment(BS18Sl_align_seqinr_ID, 
                                      matrix = "identity")
BS18Sl_tree <- nj(BS18Sl_dist)


keyvector <- as.integer(BS18Sl_tree$tip.label)


#This unshuffles the names
BS_18S_names_shuffled <- c(1:24)
for(i in 1:length(BS_18S_names)){
  BS_18S_names_shuffled[i] <- BS_18S_names[keyvector[i]]
  
}

BS18Sl_tree$tip.label <- BS_18S_names_shuffled #Renames ID labels to specimen names



plot.phylo(BS18Sl_tree, 
           type = "phylogram", 
           use.edge.length = F)


substrate_mode_og <- read.csv("Data/zo_character.csv", row.names=NULL) %>% dplyr::rename("16S" = X16S,
                                                                                   "18S" = X18S,
                                                                                   "ITS" = ITS.rDNA) 

substrate_mode <- substrate_mode_og[!is.na(substrate_mode_og$"18S"),]
#crust_mode
#substrate_mode$"Substrate"[is.na(substrate_mode$"Substrate")] <- "Unknown"
substrate_mode <- substrate_mode[!is.na(substrate_mode$"Substrate"),] 

#good_tip <- which(crust_mode$"Encrustment" != "Unknown")
substrate_mode <- substrate_mode %>% filter(Substrate != "Unknown")
#crust_mode
substrate_mode <- setNames(substrate_mode[,12],rownames(substrate_mode))

substrate_mode_shuffled <- c(1:length(substrate_mode))
for(i in 1:length(BS_18S_names)){
  substrate_mode_shuffled[i] <- substrate_mode[keyvector[i]]
  
}
cbind(substrate_mode_shuffled, BS_18S_names_shuffled)

substrate_mode <- as.factor(substrate_mode_shuffled)

plotTree(BS18Sl_tree,type="fan",fsize=0.8,ftype="i",lwd=2)
#length(levels(crust_mode))
#levels(crust_mode)
cols <- setNames(mycols,levels(substrate_mode))[c(1:length(levels(substrate_mode)))]
tiplabels(pie=to.matrix(substrate_mode,
                        levels(substrate_mode)),piecol=cols,cex=0.2)
to.matrix(substrate_mode,
          levels(substrate_mode))
add.simmap.legend(colors=cols,prompt=FALSE,x=1.5*par()$usr[3],
                  y=0*par()$usr[3],fsize=0.5)

fixed_tree <- root(BS18Sl_tree, outgroup = "Antipathozoanthus macaronesicus")
fixed_tree <- multi2di(fixed_tree)
is.rooted(fixed_tree)
is.ultrametric(fixed_tree)
is.binary.tree(fixed_tree)


branch_lengths <- branching.times(fixed_tree)
print(branch_lengths)
absolute_tree <- compute.brlen(fixed_tree, mode = "absolute") #Getting rid of all the negative branch lengths by changing them to their absolute value

fitER<-ace(substrate_mode,absolute_tree,model="ER",type="discrete")
fitER
plot.phylo(absolute_tree, type="phylogram", label.offset = 0.03, cex = 0.5)
tiplabels(pie=to.matrix(substrate_mode, levels(substrate_mode)),piecol=cols,cex=0.3)
nodelabels(node=1:absolute_tree$Nnode+Ntip(absolute_tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.3, label.offset = 3)

add.simmap.legend(colors=cols,prompt=FALSE,x=1*par()$usr[1],
                  y=3.2*par()$usr[2.5],fsize=0.8)

