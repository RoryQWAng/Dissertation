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

BS_16S <- na.omit(big_seq$"16S")

BS_18S <- na.omit(select(big_seq, "Sample", "18S"))
BS_18S_names <- BS_18S$Sample
BS_18S <- BS_18S$"18S"
BS_18S_list <- data.frame(c(1))

for(i in 1:length(BS_18S))
  BS_18S_list[[i]] <- fasta_cleaner(BS_18S[[i]], parse = F)
BS18Sl_vector <- rep(NA, length(BS_18S_list))
length(BS18Sl_vector)
for(i in 1:length(BS18Sl_vector)){
  BS18Sl_vector[i] <- BS_18S_list[[i]]
}
BS18Sl_vector_named <- Biostrings::AAStringSet(BS18Sl_vector)
 
names(BS18Sl_vector_named) <- BS_18S_names


BS18Sl_align <- msa(BS18Sl_vector_named,
                    method = "ClustalW")
BS18Sl_align_seqinr <- msaConvert(BS18Sl_align, 
                                  type = "seqinr::alignment")
BS18Sl_dist <- seqinr::dist.alignment(BS18Sl_align_seqinr, 
                                      matrix = "identity")
BS18Sl_tree <- nj(BS18Sl_dist)

plot.phylo(BS18Sl_tree, main="COI Phylogram", 
           type = "cladogram", 
           use.edge.length = T)
plot.phylo(BS18sl_tree_root, main="COI Phylogram", 
           type = "cladogram", 
           use.edge.length = T)
tiplabels()
nodelabels()
crust_mode_og <- read.csv("Data/zo_crust_2.csv", row.names=NULL) %>% dplyr::rename("16S" = X16S,
                                                                              "18S" = X18S,
                                                                              "ITS" = ITS.rDNA) 

crust_mode <- crust_mode_og[!is.na(crust_mode_og$"18S"),]
crust_mode$"Encrustment"[is.na(crust_mode$"Encrustment")] <- "Unknown"
crust_mode <- crust_mode[!is.na(crust_mode$"Encrustment"),] 
good_tip <- which(crust_mode$"Encrustment" != "Unknown")
crust_mode <- crust_mode %>% filter(Encrustment != "Unknown")
crust_mode <- setNames(crust_mode[,2],rownames(crust_mode))
crust_mode


crust_mode <- as.factor(crust_mode)
crust_mode
plotTree(BS18Sl_tree,type="fan",fsize=0.7,ftype="i",lwd=1)
length(levels(crust_mode))
levels(crust_mode)
cols <- setNames(mycols,levels(crust_mode))[c(1:length(levels(crust_mode)))]
#BS18Sl_tree$tip.label

#to.matrix(crust_mode, levels(crust_mode))
tiplabels(tip=good_tip, pie=to.matrix(crust_mode,
                        levels(crust_mode)),piecol=cols,cex=0.3)
#tiplabels(pie=to.matrix(crust_mode[BS18Sl_tree$tip.label],
                  #      levels(crust_mode)),piecol=cols,cex=0.3)
levels(crust_mode)
add.simmap.legend(colors=cols,prompt=FALSE,x=1.5*par()$usr[3],
                  y=0*par()$usr[3],fsize=0.5)
# Attempted Analysis
fitER<-ace(crust_mode,BS18Sl_tree,model="ER",type="discrete")
fitER<-ace(x,tree,model="ER",type="discrete")
fitER #Didn't like because it was not dichotomous or rooted

#######
#Redoing work but excluding NA encrustment to allow analysis
#######
#####
crust_mode2 <- crust_mode_og[!is.na(crust_mode_og$"18S"),]
crust_mode2 <- crust_mode2[!is.na(crust_mode2$"Encrustment"),]
#crust_mode2 <- na.omit(crust_mode2$Encrustment)
#crust_mode2 <- crust_mode[!is.na(crust_mode2$"Encrustment")] 
#####
#crust_mode2 <- crust_mode %>% filter(Encrustment != "Unknown")
crust_mode2 <- setNames(crust_mode2[,2],rownames(crust_mode2))
crust_mode2


crust_mode2 <- as.factor(crust_mode2)
crust_mode
plotTree(BS18Sl_tree,type="fan",fsize=0.7,ftype="i",lwd=1)
length(levels(crust_mode2))
levels(crust_mode2)
cols2 <- setNames(mycols,levels(crust_mode2))[c(1:length(levels(crust_mode2)))]
#BS18Sl_tree$tip.label



#to.matrix(crust_mode, levels(crust_mode))
tiplabels(tip=good_tip, pie=to.matrix(crust_mode2,
                                      levels(crust_mode2)),piecol=cols2,cex=0.3)
#tiplabels(pie=to.matrix(crust_mode[BS18Sl_tree$tip.label],
#      levels(crust_mode)),piecol=cols,cex=0.3)
levels(crust_mode2)
add.simmap.legend(colors=cols2,prompt=FALSE,x=1.5*par()$usr[3],
                  y=0*par()$usr[3],fsize=0.5)
###### Attempted Analysis
BS18sl_tree_root <- root(BS18Sl_tree, outgroup = "Zibrowius ammophilus")
dic_tree <- dichotomous_tree <- multi2di(BS18sl_tree_root)
BS_18S_names <- BS_18S$Sample
BS_18S <- BS_18S$"18S"
BS_18S2 <- na.omit(BS_18S$"Encrustment")
BS_18S_list <- data.frame(c(1))
#####
for(i in 1:length(BS_18S))
  BS_18S_list[[i]] <- fasta_cleaner(BS_18S[[i]], parse = F)
BS18Sl_vector <- rep(NA, length(BS_18S_list))
length(BS18Sl_vector)
for(i in 1:length(BS18Sl_vector)){
  BS18Sl_vector[i] <- BS_18S_list[[i]]
}
BS18Sl_vector_named <- Biostrings::AAStringSet(BS18Sl_vector)

names(BS18Sl_vector_named) <- BS_18S_names

fitER<-ace(crust_mode2,dic_tree,model="ER",type="discrete")

#####


BS_ITS <- na.omit(big_seq$ITS)

