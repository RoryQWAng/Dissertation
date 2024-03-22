#install.packages("bioconductor")
#install.packages("phytools")
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")

#install.packages("rentrez",dependencies = TRUE)
#install.packages("devtools")
#install.packages("ape") 
#install.packages("seqinr")
#install.packages("BiocManager")
#BiocManager::install("mxcsa")
#BiocManager::install("Biostrings")
library(devtools)
#devtools::install_github("brouwern/combio4all")
#devtools::install_github("YuLab-SMU/ggmsa")
#install_github("liamrevell/phytools")
#install.packages("tinytex")

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
detach("package:phylotools", unload=TRUE)

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


smol_seq <- read.csv("Data/small_phylo_set.csv", header=T) %>% dplyr::rename("16S" = X16S,
                                                                             "18S" = X18S,
                                                                             "ITS" = ITS.rDNA)

class(smol_seq)

smol_seq_18S <- na.omit(smol_seq$"18S")
smol_seq_18S_list <- data.frame(c(1))
for(i in 1:length(smol_seq_18S))
  smol_seq_18S_list[[i]] <- fasta_cleaner(smol_seq_18S[[i]], parse = F)

length(smol_seq_18S_list)
ss18Sl_vector <- rep(NA, length(smol_seq_18S_list))
length(ss18Sl_vector)
for(i in 1:length(ss18Sl_vector)){
  ss18Sl_vector[i] <- smol_seq_18S_list[[i]]
}


names(ss18Sl_vector) <- names(smol_seq_18S_list)  
attributes(ss18Sl_vector)
smol_seq_names <- smol_seq$Sample
smol_seq_names
names(ss18Sl_vector) <- smol_seq_names[-16]
ss18Sl_vector_ss <- Biostrings::AAStringSet(ss18Sl_vector)
names(ss18Sl_vector_ss)
smol_seq_names <- na.omit(smol_seq$"Sample")

smol_seq_loc <- na.omit(smol_seq$"Location")
smol_seq_names
smol_seq_names[-16]
names(ss18Sl_vector_ss) <- smol_seq_names[-16]

ss18Sl_align <- msa(ss18Sl_vector_ss,
                    method = "ClustalW")
names(ss18Sl_vector_ss) <- smol_seq_names[-16]
location(ss18Sl_vector_ss) <- smol_seq_loc[-16]
ss18Sl_align
typeof(ss18Sl_vector_ss)
#class(ss18Sl_align)
#is(ss18Sl_align)

class(ss18Sl_align) <- "AAMultipleAlignment"

ss18Sl_align_seqinr <- msaConvert(ss18Sl_align, 
                                  type = "seqinr::alignment")
compbio4all::print_msa(alignment = ss18Sl_align_seqinr, 
                       chunksize = 60)

class(ss18Sl_align) <- "AAMultipleAlignment"
ggmsa::ggmsa(ss18Sl_align ,  # shrooms_align, NOT shrooms_align_seqinr
             start = "180", end = "200", char_width = 1) 

ss18Sl_dist <- seqinr::dist.alignment(ss18Sl_align_seqinr, 
                                              matrix = "identity")
ss18Sl_dist_rounded <- round(ss18Sl_dist,
                                     digits = 3)
ss18Sl_tree <- nj(ss18Sl_dist_rounded)
plot.phylo(ss18Sl_tree, main="Sub Set COI Phylogram", 
           type = "phylogram", 
           use.edge.length = T)

#Using APE in an attempt to assign traits to specific samples
plot(ss18Sl_tree, type = "cladogram")
roundPhylogram(ss18Sl_tree)
plot(unroot(ss18Sl_tree),type="unrooted",no.margin=TRUE,lab4ut="axial",
     edge.width=2)

plotTree(ss18Sl_tree,offset=1)
tiplabels()
nodelabels()
plotTree(ss18Sl_tree,type="fan",fsize=0.85,lwd=1,
         ftype= "i",
         use.edge.length = F)

ss_location <- select(smol_seq[-16,], Sample) 
colour <- c("red", "green", "red", "green", "red", "green", "red", "green", "red", "green", "red", "green", "red", "green", "red", "green","red")
ss_location["colour"] = colour
ss_location <- ss_location[-16]
colnames(ss_location)
colour_vector <- unlist(colour)
factor2 <- factor(colour_vector, levels = unique(colour_vector))
combi_list <- c(ss_location, factor2)

colour.mode<-setNames(combi_list,rownames(combi_list$colour))
#feed.mode  <- feed.mode[c("colour", "Sample")]





levels(combi_list)
levels(colour.mode)
combi_list <- as.factor(combi_list)
colour.mode <- as.factor(colour.mode)
colour.mode <- read.csv("Data/fucky_wucky.csv",row.names=1)
colour.mode<- setNames(colour.mode[,1],rownames(colour.mode))
colour.mode <- as.factor(colour.mode)

plotTree(ss18Sl_tree,type="fan",fsize=0.7,ftype="i",lwd=0.2)
cols<-setNames(c("green","red"),levels(colour.mode))
tiplabels(pie=to.matrix(colour.mode[ss18Sl_tree$tip.label],
                        levels(colour.mode)),piecol=cols,cex=0.3)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

