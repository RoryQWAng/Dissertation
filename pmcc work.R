install.packages("pandas")

big_seq <- read.csv("Data/big_seq2.csv", header = T) %>% dplyr::rename("16S" = X16S,
                                                                            "18S" = X18S,
                                                                            "ITS" = ITS.rDNA)
big_seq <- select(big_seq, "Sample", "Ocean", "Canals", "Zooxanthellate", "Encrustment", "Skeleton", "Colonial", "Substrate", "Macrocnemic.Brachycnemic", "COI", "16S", "18S", "ITS") 




#For Skeleton and Zooxanthellate
big_filtered <- big_seq[complete.cases(big_seq$Zooxanthellate, big_seq$Skeleton), ]



#Changing Zooxanthellate from Y/N to 1/0
big_filtered$Zooxanthellate[big_filtered$Zooxanthellate == "Y"] <- 1
big_filtered$Zooxanthellate[big_filtered$Zooxanthellate == "N"] <- 0
#Doing the same for skeleton
big_filtered$Skeleton[na.omit(big_filtered$Skeleton == "Yes")] <- 1
big_filtered$Skeleton[na.omit(big_filtered$Skeleton == "No")] <- 0

Zooxanthellate <- as.numeric(as.character(big_filtered$Zooxanthellate))
Skeleton <- as.numeric(as.character(big_filtered$Skeleton))

cor.test(Zooxanthellate,Skeleton)


#Colonialism and Zooxanthellae, as it turns out, the only solitary zoas I couldn't find their zooxanthellism
big_filtered <- big_seq[big_seq$Colonial != "Unknown", ]
big_filtered <- big_filtered[complete.cases(big_filtered$Zooxanthellate),]

big_filtered$Zooxanthellate[big_filtered$Zooxanthellate == "Y"] <- 1
big_filtered$Zooxanthellate[big_filtered$Zooxanthellate == "N"] <- 0

big_filtered$Colonial[na.omit(big_filtered$Colonial == "Yes")] <- 1
big_filtered$Colonial[na.omit(big_filtered$Colonial == "No")] <- 0

Zooxanthellate <- as.numeric(as.character(big_filtered$Zooxanthellate))
Colonial <- as.numeric(as.character(big_filtered$Colonial))

cor.test(Zooxanthellate, Colonial)

#Colonial and Skeleton
big_filtered <- big_seq[big_seq$Colonial != "Unknown", ]
big_filtered <- big_filtered[complete.cases(big_filtered$Skeleton),]

big_filtered$Skeleton[big_filtered$Skeleton == "Yes"] <- 1
big_filtered$Skeleton[big_filtered$Skeleton == "No"] <- 0

big_filtered$Colonial[na.omit(big_filtered$Colonial == "Yes")] <- 1
big_filtered$Colonial[na.omit(big_filtered$Colonial == "No")] <- 0

Skeleton <- as.numeric(as.character(big_filtered$Skeleton))
Colonial <- as.numeric(as.character(big_filtered$Colonial))

cor.test(Skeleton, Colonial)

#Colonial and M/B
big_filtered <- big_seq[big_seq$Macrocnemic.Brachycnemic != "Unkown", ]
big_filtered <- big_filtered[big_filtered$Macrocnemic.Brachycnemic != "Unknown", ]
big_filtered <- big_filtered[complete.cases(big_filtered$Colonial),]

big_filtered$Macrocnemic.Brachycnemic[big_filtered$Macrocnemic.Brachycnemic == "Macrocnemic"] <- 1
big_filtered$Macrocnemic.Brachycnemic[big_filtered$Macrocnemic.Brachycnemic == "Brachycnemic"] <- 0

big_filtered$Skeleton[na.omit(big_filtered$Skeleton == "Yes")] <- 1
big_filtered$Skeleton[na.omit(big_filtered$Skeleton == "No")] <- 0

Skeleton <- as.numeric(as.character(big_filtered$Skeleton))
MB <- as.numeric(as.character(big_filtered$Macrocnemic.Brachycnemic))

cor.test(MB, Skeleton)

#Work on Encurstment - Encrustment and Skeleton
big_filtered <- big_seq[big_seq$Encrustment != "Yes", ]
big_filtered <- big_filtered[complete.cases(big_filtered$Encrustment),]
big_filtered <- big_filtered[complete.cases(big_filtered$Skeleton),]

big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "None")] <- 0
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Endoderm")] <- 1
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Centre Mesoglea")] <- 2
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Outer Mesoglea")] <- 3
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Ectoderm")] <- 4
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Ectoderm/Mesoglea")] <- 5
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Ectoderm")] <- 6
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Polyp")] <- 7

big_filtered$Skeleton[na.omit(big_filtered$Skeleton == "Yes")] <- 1
big_filtered$Skeleton[na.omit(big_filtered$Skeleton == "No")] <- 0

Skeleton <- as.numeric(as.character(big_filtered$Skeleton))
Crust <- as.numeric(as.character(big_filtered$Encrustment))

cor.test(Skeleton, Crust)

#BM and Encrustment
big_filtered <- big_seq[complete.cases(big_filtered$Encrustment),]
big_filtered <- big_filtered[complete.cases(big_filtered$Macrocnemic.Brachycnemic),]
big_filtered <- big_filtered[complete.cases(big_filtered$Encrustment),]
big_filtered <- big_filtered[big_filtered$Encrustment != "Yes", ]
big_filtered <- big_filtered[big_filtered$Macrocnemic.Brachycnemic != "Unknown", ]
big_filtered <- big_seq[big_seq$Macrocnemic.Brachycnemic != "Unkown", ]

big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "None")] <- 0
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Endoderm")] <- 1
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Centre Mesoglea")] <- 2
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Outer Mesoglea")] <- 3
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Ectoderm")] <- 4
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Ectoderm/Mesoglea")] <- 5
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Ectoderm")] <- 6
big_filtered$Encrustment[na.omit(big_filtered$Encrustment == "Polyp")] <- 7

big_filtered$Macrocnemic.Brachycnemic[big_filtered$Macrocnemic.Brachycnemic == "Macrocnemic"] <- 1
big_filtered$Macrocnemic.Brachycnemic[big_filtered$Macrocnemic.Brachycnemic == "Brachycnemic"] <- 0

MB <- as.numeric(as.character(big_filtered$Macrocnemic.Brachycnemic))
Crust <- as.numeric(as.character(big_filtered$Encrustment))
