#Préparation de données et comparaison de mocks 16S
#Mathias

####Import OTU/TAX####

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

####Renommage####

#Renommage des mocks
colnames(otu_16S_vsearch)[colnames(otu_16S_vsearch)=="Mock"] <- "Mock_vsearch"
colnames(otu_16S_deblur)[colnames(otu_16S_deblur)=="Mock"] <- "Mock_deblur"
colnames(otu_16S_dada2)[colnames(otu_16S_dada2)=="Mock"] <- "Mock_dada2"

#Renommage des Taxon
names(taxa_16S_deblur)[1]<-"OTUID"
colnames(taxa_16S_vsearch)[colnames(taxa_16S_vsearch)=="Taxon"] <- "Taxon_vsearch"
colnames(taxa_16S_deblur)[colnames(taxa_16S_deblur)=="Taxon"] <- "Taxon_deblur"
colnames(taxa_16S_dada2)[colnames(taxa_16S_dada2)=="Taxon"] <- "Taxon_dada2"

####Import SEQ####

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

####Merge####

#Merge des tables pour avoir comme colonnes finales OTU/Seq/Mocks
seq_merge_16S_vsearch<-merge(s_vsearch[,c("OTUID","sequence_vsearch")],otu_16S_vsearch[,c("OTUID","Mock_vsearch")],by.x = c("OTUID"),all.x=F, all.y=F)
seq_merge_16S_deblur<-merge(s_deblur[,c("OTUID","sequence_deblur")],otu_16S_deblur[,c("OTUID","Mock_deblur")],by.x = c("OTUID"),all.x=F, all.y=F)
seq_merge_16S_dada2<-merge(s_dada2[,c("OTUID","sequence_dada2")],otu_16S_dada2[,c("OTUID","Mock_dada2")],by.x = c("OTUID"),all.x=F, all.y=F)

seq_merge_vsearch_deblur<-merge(seq_merge_16S_vsearch[,c("sequence_vsearch","Mock_vsearch")],seq_merge_16S_deblur[,c("sequence_deblur","Mock_deblur")],by.x = c("sequence_vsearch"),by.y = c("sequence_deblur"),all.x=T, all.y=T)
seq_merge_all<-merge(seq_merge_vsearch_deblur[,c("sequence_vsearch","Mock_vsearch","Mock_deblur")],seq_merge_16S_dada2[,c("sequence_dada2","Mock_dada2")],by.x = c("sequence_vsearch"),by.y = c("sequence_dada2"),all.x=T, all.y=T)

####Manip####

#Renommage de la colonne des séquences
colnames(seq_merge_all)[colnames(seq_merge_all)=="sequence_vsearch"] <- "sequence"

#Remplacement des NA avec des 0
seq_merge_all[is.na(seq_merge_all)] <- 0

#Supression des lignes où il y a que des comptes d'OTU à 0
seq_final <- seq_merge_all[!rowSums(seq_merge_all[, -1] == 0) == (ncol(seq_merge_all)-1), ]

#Merge des OTUIDs pour récupérer ensuite la taxonomie
#seq_final<-merge(seq_final[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2")],seq_merge_16S_vsearch[,c("OTUID","sequence_vsearch")],by.x = c("sequence"),by.y = c("sequence_vsearch"),all.x=T, all.y=F)
#seq_final<-merge(seq_final[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2","OTUID")],seq_merge_16S_deblur[,c("OTUID","sequence_deblur")],by.x = c("sequence"),by.y = c("sequence_deblur"),all.x=T, all.y=F)

tutu<-seq_merge_16S_vsearch
tata<-seq_merge_16S_deblur
colnames(tutu)[colnames(tutu)=="sequence_vsearch"] <- "sequence"
colnames(tutu)[colnames(tutu)=="Mock_vsearch"] <- "Mock"
colnames(tata)[colnames(tata)=="sequence_deblur"] <- "sequence"
colnames(tata)[colnames(tata)=="Mock_deblur"] <- "Mock"
test<-rbind(tutu,tata)
final<-merge(seq_final[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2")],test[,c("OTUID","sequence")],by.x = c("sequence"),all.x=T, all.y=F)

ruru<-taxa_16S_deblur
roro<-taxa_16S_vsearch
colnames(ruru)[colnames(ruru)=="Taxon_deblur"] <- "Taxon"
colnames(roro)[colnames(roro)=="Taxon_vsearch"] <- "Taxon"
final2<-rbind(ruru,roro)
final2<-merge(final[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2")],final2[,c("OTUID","Taxon")],by.x = c("OTUID"),all.x=T, all.y=F)


df2 <- data.frame(newCol=paste(final2$OTUID,final2$Taxon,sep=";"))
names(df2)[1]<-"OTUID"

#Recherche pour séparer la colonnes Taxon en différentes colonnes pour ressembler à la nommenclature du tuto
#"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
library(stringr)
test <- str_split_fixed(df2$OTUID,";D_.__",8)
test<-as.data.frame(test)
names(test) <- c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
colnames(test)

fin<-merge(test[,c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")],final2[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2")],by.x = c("OTUID"),all.x=T, all.y=F)
fin<-fin[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
fin<-rownames("OTUID")


#Export in .tsv file
write.table(x = fin, file = "/home/galati/Téléchargements/mock_table_16S.tsv",sep="\t",dec=",",row.names = F)
write.table(x = seq_final, file = "/home/galati/Téléchargements/seq_final.tsv",sep="\t",dec=",",row.names = F)

#Visualisation des données attendues
library(phyloseq)
data(GlobalPatterns)
colnames(tax_table(GlobalPatterns))


"______________________________________________________________________________________________________________________________________________________"



###Anciennes manips

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


"______________________________________________________________________________________________________________________________________________________"


otu_seqs <- colnames(fin)
otu_headers <- vector(dim(fin)[2], mode="character")

for (i in 1:dim(fin)[2]) {
  otu_headers[i] <- paste(">otu", i, sep="_")
}

# making and writing out a fasta of our final otu seqs:
otu_fasta <- c(rbind(otu_headers, otu_seqs))
#write(otu_fasta, "otus.fa")

# count table:
otu_tab <- t(fin)
row.names(otu_tab) <- sub(">", "", otu_headers)

# tax table:
otu_tax <- taxa
row.names(otu_tax) <- sub(">", "", otu_headers)
