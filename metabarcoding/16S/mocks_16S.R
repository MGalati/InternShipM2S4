#Préparation de données et comparaison de mocks 16S
#Mathias

#########

# Importation 16S vsearch
path = "/home/galati/Téléchargements/export_16S_vsearch/export"
setwd(path)
otu_16S_vsearch <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE)
taxa_16S_vsearch <- read.table(file = '16S_vsearch_taxonomy.tsv', sep = '\t', dec=".", header = TRUE)

# Importation 16S deblur
path = "/home/galati/Téléchargements/export_16S_deblur/export"
setwd(path)
otu_16S_deblur <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE)
taxa_16S_deblur <- read.table(file = '16S_deblur_taxonomy.tsv', sep = '\t', dec=".", header = TRUE)

# Importation 16S dada2
path = "/home/galati/Téléchargements/export_16S_dada2/export"
setwd(path)
otu_16S_dada2 <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE)
taxa_16S_dada2 <- read.table(file = '16S_dada2_taxonomy.tsv', sep = '\t', dec=".", header = TRUE)

#########

#Renommage des mocks
colnames(otu_16S_vsearch)[colnames(otu_16S_vsearch)=="Mock"] <- "Mock_vsearch"
colnames(otu_16S_deblur)[colnames(otu_16S_deblur)=="Mock"] <- "Mock_deblur"
colnames(otu_16S_dada2)[colnames(otu_16S_dada2)=="Mock"] <- "Mock_dada2"

#Renommage des Taxon
colnames(taxa_16S_vsearch)[colnames(taxa_16S_vsearch)=="Taxon"] <- "Taxon_vsearch"
colnames(taxa_16S_deblur)[colnames(taxa_16S_deblur)=="Taxon"] <- "Taxon_deblur"
colnames(taxa_16S_dada2)[colnames(taxa_16S_dada2)=="Taxon"] <- "Taxon_dada2"

#########

library("Biostrings")

#Importation des séquences en dataframe pour les OTUID de vsearch
path = "/home/galati/Téléchargements/export_16S_vsearch/export"
setwd(path)
fastaFile <- readDNAStringSet("dna-sequences.fasta")
OTUID = names(fastaFile)
sequence_vsearch = paste(fastaFile)
s_vsearch <- data.frame(OTUID, sequence_vsearch)

#Importation des séquences en dataframe pour les OTUID de deblur
path = "/home/galati/Téléchargements/export_16S_deblur/export"
setwd(path)
fastaFile <- readDNAStringSet("dna-sequences.fasta")
OTUID = names(fastaFile)
sequence_deblur = paste(fastaFile)
s_deblur <- data.frame(OTUID, sequence_deblur)

#Importation des séquences en dataframe pour les OTUID de dada2
path = "/home/galati/Téléchargements/export_16S_dada2/export"
setwd(path)
fastaFile <- readDNAStringSet("dna-sequences.fasta")
OTUID = names(fastaFile)
sequence_dada2 = paste(fastaFile)
s_dada2 <- data.frame(OTUID, sequence_dada2)

#########

#Merge des tables pour avoir comme colonnes finales OTU/Seq/Mocks
seq_merge_16S_vsearch<-merge(s_vsearch[,c("OTUID","sequence_vsearch")],otu_16S_vsearch[,c("OTUID","Mock_vsearch")],by.x = c("OTUID"),all.x=F, all.y=F)
seq_merge_16S_deblur<-merge(s_deblur[,c("OTUID","sequence_deblur")],otu_16S_deblur[,c("OTUID","Mock_deblur")],by.x = c("OTUID"),all.x=F, all.y=F)
seq_merge_16S_dada2<-merge(s_dada2[,c("OTUID","sequence_dada2")],otu_16S_dada2[,c("OTUID","Mock_dada2")],by.x = c("OTUID"),all.x=F, all.y=F)

seq_merge_vsearch_deblur<-merge(seq_merge_16S_vsearch[,c("sequence_vsearch","Mock_vsearch")],seq_merge_16S_deblur[,c("sequence_deblur","Mock_deblur")],by.x = c("sequence_vsearch"),by.y = c("sequence_deblur"),all.x=T, all.y=T)
seq_merge_all<-merge(seq_merge_vsearch_deblur[,c("sequence_vsearch","Mock_vsearch","Mock_deblur")],seq_merge_16S_dada2[,c("sequence_dada2","Mock_dada2")],by.x = c("sequence_vsearch"),by.y = c("sequence_dada2"),all.x=T, all.y=T)

#########

#Renommage de la colonne des séquences
colnames(seq_merge_all)[colnames(seq_merge_all)=="sequence_vsearch"] <- "sequence"

#Remplacement des NA avec des 0
seq_merge_all[is.na(seq_merge_all)] <- 0

#Supression des lignes où il y a que des comptes d'OTU à 0
seq_final <- seq_merge_all[!rowSums(seq_merge_all[, -1] == 0) == (ncol(seq_merge_all)-1), ]

#Export .csv file
write.table(x = seq_final, file = "/home/galati/Téléchargements/mock_table_16S.tsv")



"______________________________________________________________________________________________________________________________________________________"



###Anciennes manips

# Merge fichiers vsearch
OTU_merge_16S_vsearch<-merge(taxa_16S_vsearch[,c("OTUID","Taxon_vsearch")],otu_16S_vsearch[,c("OTUID","Mock_vsearch")],by.x = c("OTUID"),all.x=F, all.y=F)

# Merge fichiers deblur
OTU_merge_16S_deblur<-merge(taxa_16S_deblur[,c("OTUID","Taxon_deblur")],otu_16S_deblur[,c("OTUID","Mock_deblur")],by.x = c("OTUID"),all.x=F, all.y=F)

# Merge fichiers dada2
OTU_merge_16S_dada2<-merge(taxa_16S_dada2[,c("OTUID","Taxon_dada2")],otu_16S_dada2[,c("OTUID","Mock_dada2")],by.x = c("OTUID"),all.x=F, all.y=F)

# Merge vsearch/deblur
OTU_merge_vsearch_deblur<-merge(OTU_merge_16S_vsearch[,c("OTUID","Taxon_vsearch","Mock_vsearch")],OTU_merge_16S_deblur[,c("OTUID","Taxon_deblur","Mock_deblur")],by = c("OTUID"),all.x=T, all.y=T)

# Merge vsearch/deblur/dada2
table <-merge(OTU_merge_vsearch_deblur[,c("OTUID","Taxon_vsearch","Taxon_deblur","Mock_vsearch","Mock_deblur")],OTU_merge_16S_dada2[,c("OTUID","Taxon_dada2","Mock_dada2")],by = c("OTUID"),all.x=T, all.y=T)

#Réarrangement des colonnes
table <- table[,c("OTUID","Taxon_vsearch","Taxon_deblur","Taxon_dada2","Mock_vsearch","Mock_deblur","Mock_dada2")]

#OTUID en première colonne
#library(textshape)
#table <- column_to_rownames(table, 'OTUID')

#Remplacement des NA avec des 0
table[is.na(table)] <- 0

#Suppression des lignes ayant que des comptes de 0 séquence OTU détecté quand il n'y a que les OTUID en colonne 1 puis les Mocks
#table1 <- table[!rowSums(table[, -1] == 0) == (ncol(table)-1), ]

#Crée une table avec OTUID et les 3 Mocks
tableT <- table[,c(1,5:7)]
tableT <- tableT[!rowSums(tableT == 0) == (ncol(tableT)-1), ]

#Concaténation des données
final <-merge(tableT[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2")],table[,c("OTUID","Taxon_vsearch","Taxon_deblur","Taxon_dada2")],by.x = c("OTUID"),all.x=T, all.y=F)



"______________________________________________________________________________________________________________________________________________________"



#Merge des séquences aux OTUID 
merge_seq_vsearch<-merge(final[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","Taxon_vsearch","Taxon_deblur","Taxon_dada2")],s_vsearch[,c("OTUID","sequence_vsearch")],by.x = c("OTUID"),all.x=T, all.y=F)
merge_seq_vsearch_deblur<-merge(merge_seq_vsearch[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","Taxon_vsearch","Taxon_deblur","Taxon_dada2","sequence_vsearch")],s_deblur[,c("OTUID","sequence_deblur")],by.x = c("OTUID"),all.x=T, all.y=F)
merge_seq_all<-merge(merge_seq_vsearch_deblur[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","Taxon_vsearch","Taxon_deblur","Taxon_dada2","sequence_vsearch","sequence_deblur")],s_dada2[,c("OTUID","sequence_dada2")],by.x = c("OTUID"),all.x=T, all.y=F)


"______________________________________________________________________________________________________________________________________________________"





