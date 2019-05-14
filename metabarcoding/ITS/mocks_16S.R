#Comparaison de mocks
#Mathias


"______________________________________________________________________________________________________________________________________________________"


# Importation 16S vsearch
path = "/home/galati/Téléchargements/export_16S_vsearch/export"
setwd(path)
otu_16S_vsearch <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE, stringsAsFactors = FALSE)
taxa_16S_vsearch <- read.table(file = '16S_vsearch_taxonomy.tsv', sep = '\t', dec=".", header = TRUE, stringsAsFactors = FALSE)

# Importation 16S deblur
path = "/home/galati/Téléchargements/export_16S_deblur/export"
setwd(path)
otu_16S_deblur <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE, stringsAsFactors = FALSE)
taxa_16S_deblur <- read.table(file = '16S_deblur_taxonomy.tsv', sep = '\t', dec=".", header = TRUE, stringsAsFactors = FALSE)

# Importation 16S dada2
path = "/home/galati/Téléchargements/export_16S_dada2/export"
setwd(path)
otu_16S_dada2 <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE, stringsAsFactors = FALSE)
taxa_16S_dada2 <- read.table(file = '16S_dada2_taxonomy.tsv', sep = '\t', dec=".", header = TRUE, stringsAsFactors = FALSE)


"______________________________________________________________________________________________________________________________________________________"


#Renommage des mocks
colnames(otu_16S_vsearch)[colnames(otu_16S_vsearch)=="Mock"] <- "Mock_vsearch"
colnames(otu_16S_deblur)[colnames(otu_16S_deblur)=="Mock"] <- "Mock_deblur"
colnames(otu_16S_dada2)[colnames(otu_16S_dada2)=="Mock"] <- "Mock_dada2"

#Renommage des Taxon
colnames(taxa_16S_vsearch)[colnames(taxa_16S_vsearch)=="Taxon"] <- "Taxon_vsearch"
colnames(taxa_16S_deblur)[colnames(taxa_16S_deblur)=="Taxon"] <- "Taxon_deblur"
colnames(taxa_16S_dada2)[colnames(taxa_16S_dada2)=="Taxon"] <- "Taxon_dada2"


"______________________________________________________________________________________________________________________________________________________"


# Merge fichiers vsearch
OTU_merge_16S_vsearch<-merge(taxa_16S_vsearch[,c("OTUID","Taxon_vsearch")],otu_16S_vsearch[,c("OTUID","Mock_vsearch")],by.x = c("OTUID"),all.x=F, all.y=F)

# Merge fichiers deblur
OTU_merge_16S_deblur<-merge(taxa_16S_deblur[,c("OTUID","Taxon_deblur")],otu_16S_deblur[,c("OTUID","Mock_deblur")],by.x = c("OTUID"),all.x=F, all.y=F)

# Merge fichiers dada2
OTU_merge_16S_dada2<-merge(taxa_16S_dada2[,c("OTUID","Taxon_dada2")],otu_16S_dada2[,c("OTUID","Mock_dada2")],by.x = c("OTUID"),all.x=F, all.y=F)


"______________________________________________________________________________________________________________________________________________________"


# Merge vsearch/deblur
OTU_merge_vsearch_deblur<-merge(OTU_merge_16S_vsearch[,c("OTUID","Taxon_vsearch","Mock_vsearch")],OTU_merge_16S_deblur[,c("OTUID","Taxon_deblur","Mock_deblur")],by = c("OTUID"),all.x=T, all.y=T)

# Merge vsearch/deblur/dada2
table <-merge(OTU_merge_vsearch_deblur[,c("OTUID","Taxon_vsearch","Taxon_deblur","Mock_vsearch","Mock_deblur")],OTU_merge_16S_dada2[,c("OTUID","Taxon_dada2","Mock_dada2")],by = c("OTUID"),all.x=T, all.y=T)


"______________________________________________________________________________________________________________________________________________________"


#Réarrangement des colonnes
table <- table[,c("OTUID","Taxon_vsearch","Taxon_deblur","Taxon_dada2","Mock_vsearch","Mock_deblur","Mock_dada2")]

#OTUID en première colonne
#library(textshape)
#table <- column_to_rownames(table, 'OTUID')

#Remplacement des NA avec des 0
table[is.na(table)] <- 0

#Suppression des lignes ayant que des comptes de 0 séquence OTU détecté quand il n'y a que les OTUID en colonne 1 puis les Mocks
#table1 <- table[!rowSums(table[, -1] == 0) == (ncol(table)-1), ]

#Work in progress = Supression des lignes ayant pour les colonnes Mocks une somme de 0+0+0.
#Essai 1
table2 <- as.matrix(table)
table2[,5:7] <- as.integer(table2[,5:7])
for (i in 5:ncol(table2)){
  if (sum(table2[,i])==0){
    print(names(table[i]))
  }
}

#Essai 2
tableT <- table[!rowSums(table[, -1:-4] == 0) == (ncol(table)-1:-4), ]


"______________________________________________________________________________________________________________________________________________________"


#Export .csv file
write.table(x = table, file = "/home/galati/Téléchargements/mock_table_16S.tsv")
