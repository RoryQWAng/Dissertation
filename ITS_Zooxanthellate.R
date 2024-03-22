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
big_seq <- read.csv("Data/Zoa_character.csv", header = T) %>% dplyr::rename("16S" = X16S,
                                                                            "18S" = X18S,
                                                                            "ITS" = ITS.rDNA) 
big_seq <- select(big_seq, "Sample", "Canals", "Zooxanthellate", "Gender", "Encrustment","Skeleton", "COI", "16S", "18S", "ITS")
BS_ITS <- na.omit(select(big_seq, "Sample", "ITS", "Zooxanthellate"))
BS_ITS <- select(BS_ITS, "Sample", "ITS") 
BS_ITS_names <- BS_ITS$"Sample"
BS_ITS$"ID" <- 1:nrow(BS_ITS)

BS_ITS_ID <- BS_ITS$ID
BS_ITS <- BS_ITS$"ITS"
BS_ITS_list <- data.frame(c(1))

for(i in 1:length(BS_ITS))
  BS_ITS_list[[i]] <- fasta_cleaner(BS_ITS[[i]], parse = F)
BSITSl_vector <- rep(NA, length(BS_ITS_list))
length(BSITSl_vector)
for(i in 1:length(BSITSl_vector)){
  BSITSl_vector[i] <- BS_ITS_list[[i]]
}

BSITSl_vector_ID <- Biostrings::AAStringSet(BSITSl_vector)


names(BSITSl_vector_ID) <- BS_ITS_ID



BSITSl_align_ID <- msa(BSITSl_vector_ID,
                       method = "ClustalW")

BSITSl_align_seqinr_ID <- msaConvert(BSITSl_align_ID, 
                                     type = "seqinr::alignment")

BSITSl_dist <- seqinr::dist.alignment(BSITSl_align_seqinr_ID, 
                                      matrix = "identity")
BSITSl_tree <- nj(BSITSl_dist)


keyvector <- as.integer(BSITSl_tree$tip.label)


#This unshuffles the names
BS_ITS_names_shuffled <- c(1:24)
for(i in 1:length(BS_ITS_names)){
  BS_ITS_names_shuffled[i] <- BS_ITS_names[keyvector[i]]
  
}

BSITSl_tree$tip.label <- BS_ITS_names_shuffled #Renames ID labels to specimen names



plot.phylo(BSITSl_tree, 
           type = "phylogram", 
           use.edge.length = F)


crust_mode_og <- read.csv("Data/zo_character.csv", row.names=NULL) %>% dplyr::rename("16S" = X16S,
                                                                                   "18S" = X18S,
                                                                                   "ITS" = ITS.rDNA) 

crust_mode <- crust_mode_og[!is.na(crust_mode_og$"ITS"),]
#crust_mode
#crust_mode$"Zooxanthellate"[is.na(crust_mode$"Zooxanthellate")] <- "Unknown"
crust_mode <- crust_mode[!is.na(crust_mode$"Zooxanthellate"),] 

#good_tip <- which(crust_mode$"Zooxanthellate" != "Unknown")
#crust_mode <- crust_mode %>% filter(Zooxanthellate != "Unknown")
#crust_mode
crust_mode <- setNames(crust_mode[,7],rownames(crust_mode))

crust_mode_shuffled <- c(1:length(BS_ITS_names))
for(i in 1:length(BS_ITS_names)){
  crust_mode_shuffled[i] <- crust_mode[keyvector[i]]
  
}
cbind(crust_mode_shuffled, BS_ITS_names_shuffled)

crust_mode <- as.factor(crust_mode_shuffled)

plotTree(BSITSl_tree,type="phylogram",fsize=0.8,ftype="i",lwd=2)
#length(levels(crust_mode))
#levels(crust_mode)
cols <- setNames(mycols,levels(crust_mode))[c(1:length(levels(crust_mode)))]
tiplabels(pie=to.matrix(crust_mode,
                        levels(crust_mode)),piecol=cols,cex=0.2)
to.matrix(crust_mode,
          levels(crust_mode))
add.simmap.legend(colors=cols,prompt=FALSE,x=1.5*par()$usr[3],
                  y=0*par()$usr[3],fsize=0.5)

fixed_tree <- root(BSITSl_tree, outgroup = "Parazoanthus sp. New Caledonia Shallow")
fixed_tree <- multi2di(fixed_tree)
is.rooted(fixed_tree)
is.ultrametric(fixed_tree)
is.binary.tree(fixed_tree)


branch_lengths <- branching.times(fixed_tree)
print(branch_lengths)
absolute_tree <- compute.brlen(fixed_tree, mode = "absolute") #Getting rid of all the negative branch lengths by changing them to their absolute value

fitER<-ace(crust_mode,absolute_tree,model="ER",type="discrete")
fitER
plot.phylo(absolute_tree, type="phylogram", label.offset = 0.03, cex = 0.5)
tiplabels(pie=to.matrix(crust_mode, levels(crust_mode)),piecol=cols,cex=0.3)
nodelabels(node=1:absolute_tree$Nnode+Ntip(absolute_tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.3, label.offset = 3)

add.simmap.legend(colors=cols,prompt=FALSE,x=1*par()$usr[1],
                  y=3.2*par()$usr[2.5],fsize=0.8)

