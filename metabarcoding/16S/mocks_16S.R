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

####Récupération des séquences pour avoir une table combinant les séquences et les identifiants Qiime2####
colnames(seq_merge_16S_vsearch)[colnames(seq_merge_16S_vsearch)=="sequence_vsearch"] <- "sequence"
colnames(seq_merge_16S_vsearch)[colnames(seq_merge_16S_vsearch)=="Mock_vsearch"] <- "Mock"
colnames(seq_merge_16S_deblur)[colnames(seq_merge_16S_deblur)=="sequence_deblur"] <- "sequence"
colnames(seq_merge_16S_deblur)[colnames(seq_merge_16S_deblur)=="Mock_deblur"] <- "Mock"
otu_vsearch_deblur<-rbind(seq_merge_16S_vsearch,seq_merge_16S_deblur)
otu_seq_id<-merge(seq_final[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2")],otu_vsearch_deblur[,c("OTUID","sequence")],by.x = c("sequence"),all.x=T, all.y=F)

####Obtention d'une table OTU/Mocks/Taxonomie####
colnames(taxa_16S_deblur)[colnames(taxa_16S_deblur)=="Taxon_deblur"] <- "Taxon"
colnames(taxa_16S_vsearch)[colnames(taxa_16S_vsearch)=="Taxon_vsearch"] <- "Taxon"
colnames(taxa_16S_vsearch)[colnames(taxa_16S_vsearch)=="Feature.ID"] <- "OTUID"
colnames(taxa_16S_dada2)[colnames(taxa_16S_dada2)=="Taxon_dada2"] <- "Taxon"
colnames(taxa_16S_dada2)[colnames(taxa_16S_dada2)=="Feature.ID"] <- "OTUID"
tax_tot<-rbind(taxa_16S_deblur,taxa_16S_vsearch, taxa_16S_dada2)
otu_mock_tax<-merge(otu_seq_id[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2")],tax_tot[,c("OTUID","Taxon")],by.x = c("OTUID"),all.x=T, all.y=T)

####Combiner les infos d'OTU et de Taxonomie####
comb_tax <- data.frame(newCol=paste(otu_mock_tax$OTUID,otu_mock_tax$Taxon,sep=";"))
names(comb_tax)[1]<-"OTUID"

####Recherche pour séparer la colonnes Taxon en différentes colonnes pour ressembler à la nommenclature du tuto####
#"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
library(stringr)
tax_split <- str_split_fixed(comb_tax$OTUID,";D_.__",8)
tax_split<-as.data.frame(tax_split)
names(tax_split) <- c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
#Vérification
colnames(tax_split)

id_mocks_taxsplit<-merge(tax_split[,c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")],otu_mock_tax[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2")],by.x = c("OTUID"),all.x=T, all.y=F)
id_mocks_taxsplit<-id_mocks_taxsplit[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","Kingdom","Phylum","Class","Order","Family","Genus","Species")]


####True Mock#### 
id_mocks_taxsplit$True_OTU<-0
id_mocks_taxsplit<-id_mocks_taxsplit[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","True_OTU","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
id_mocks_taxsplit$OTUID<-as.character(id_mocks_taxsplit$OTUID)
id_mocks_taxsplit$Order<-as.character(id_mocks_taxsplit$Order)
id_mocks_taxsplit$Family<-as.character(id_mocks_taxsplit$Family)
id_mocks_taxsplit$Genus<-as.character(id_mocks_taxsplit$Genus)
id_mocks_taxsplit$Species<-as.character(id_mocks_taxsplit$Species)

mock1<-c("Acinetobacter_baumannii",0,0,0,189,"Bacteria","Proteobacteria","Gammaproteobacteria","Pseudomonadales","Moraxellaceae","Acinetobacter","Acinetobacter_baumannii")
mock2<-c("Actinomyces_odontolyticus",0,0,0,189,"Bacteria","Actinobacteria","Actinobacteria","Actinomycetales","Actinomycetaceae","Actinomyces","Actinomyces_odontolyticus")
mock3<-c("Bacillus_cereus",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Bacillales","Bacillaceae","Bacillus","Bacillus_cereus")
mock4<-c("Bacteroides_vulgatus",0,0,0,189,"Bacteria","Bacteroidetes","Bacteroidia","Bacteroidales","Bacteroidaceae","Bacteroides","Bacteroides_vulgatus")
mock5<-c("Clostridium_beijerinckii",0,0,0,189,"Bacteria","Firmicutes","Clostridia","Clostridiales","Clostridiaceae","Clostridium","Clostridium_beijerinckii")
mock6<-c("Deinococcus_radiodurans",0,0,0,189,"Bacteria","Deinococcus-Thermus","Deinococci","Deinococcales","Deinococcaceae","Deinococcus","Deinococcus_radiodurans")
mock7<-c("Enterococcus_faecalis",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Lactobacillales","Enterococcaceae","Enterococcus","Enterococcus_faecalis")
mock8<-c("Escherichia_coli",0,0,0,189,"Bacteria","Proteobacteria","Gammaproteobacteria","Enterobacteriales","Enterobacteriaceae","Escherichia","Escherichia_coli")
mock9<-c("Helicobacter_pylori",0,0,0,189,"Bacteria","Epsilonbacteraeota","Campylobacteria","Campylobacterales","Helicobacteraceae","Helicobacter","Helicobacter_pylori")
mock10<-c("Lactobacillus_gasseri",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Lactobacillales","Lactobacillaceae","Lactobacillus","Lactobacillus_gasseri")
mock11<-c("Listeria_monocytogenes",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Bacillales","Listeriaceae","Listeria","Listeria_monocytogenes")
mock12<-c("Neisseria_meningitidis",0,0,0,189,"Bacteria","Proteobacteria","Gammaproteobacteria","Betaproteobacteriales","Neisseriaceae","Neisseria","Neisseria_meningitidis")
mock13<-c("Porphyromonas_gingivalis",0,0,0,189,"Bacteria","Bacteroidetes","Bacteroidia","Bacteroidales","Bacteroidaceae","Bacteroides","Porphyromonas_gingivalis")
mock14<-c("Propionibacterium_acnes",0,0,0,189,"Bacteria","Actinobacteria","Actinobacteria","Propionibacteriales","Propionibacteriaceae","Cutibacterium","Propionibacterium_acnes")
mock15<-c("Pseudomonas_aeruginosa",0,0,0,189,"Bacteria","Proteobacteria","Gammaproteobacteria","Pseudomonadales","Pseudomonadaceae","Pseudomonas","Pseudomonas_aeruginosa")
mock16<-c("Rhodobacter_sphaeroides",0,0,0,189,"Bacteria","Proteobacteria","Alphaproteobacteria","Rhodobacterales","Rhodobacteraceae","Rhodobacter","Rhodobacter_sphaeroides")
mock17<-c("Staphylococcus_aureus",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Bacillales","Staphylococcaceae","Staphylococcus","Staphylococcus_aureus")
mock18<-c("Staphylococcus_epidermidis",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Bacillales","Staphylococcaceae","Staphylococcus","Staphylococcus_epidermidis")
mock19<-c("Streptococcus_agalactiae",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Lactobacillales","Streptococcaceae","Streptococcus","Streptococcus_agalactiae")
mock20<-c("Streptococcus_mutans",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Lactobacillales","Streptococcaceae","Streptococcus","Streptococcus_mutans")
mock21<-c("Streptococcus_pneumoniae",0,0,0,189,"Bacteria","Firmicutes","Bacilli","Lactobacillales","Streptococcaceae","Streptococcus","Streptococcus_pneumoniae")
id_mocks_taxsplit<-rbind(id_mocks_taxsplit,mock1,mock2,mock3,mock4,mock5,mock6,mock7,mock8,mock9,mock10,mock11,mock12,mock13,mock14,mock15,mock16,mock17,mock18,mock19,mock20,mock21)

####Création des tables pour phyloseq####
#Table OTU
otu<-otu_seq_id[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2")]
otu$sequence<-as.character(otu$sequence)
otu$True_OTU<-0
tax1<-c("Acinetobacter_baumannii",0,0,0,189)
tax2<-c("Actinomyces_odontolyticus",0,0,0,189)
tax3<-c("Bacillus_cereus",0,0,0,189)
tax4<-c("Bacteroides_vulgatus",0,0,0,189)
tax5<-c("Clostridium_beijerinckii",0,0,0,189)
tax6<-c("Deinococcus_radiodurans",0,0,0,189)
tax7<-c("Enterococcus_faecalis",0,0,0,189)
tax8<-c("Escherichia_coli",0,0,0,189)
tax9<-c("Helicobacter_pylori",0,0,0,189)
tax10<-c("Lactobacillus_gasseri",0,0,0,189)
tax11<-c("Listeria_monocytogenes",0,0,0,189)
tax12<-c("Neisseria_meningitidis",0,0,0,189)
tax13<-c("Porphyromonas_gingivalis",0,0,0,189)
tax14<-c("Propionibacterium_acnes",0,0,0,189)
tax15<-c("Pseudomonas_aeruginosa",0,0,0,189)
tax16<-c("Rhodobacter_sphaeroides",0,0,0,189)
tax17<-c("Staphylococcus_aureus",0,0,0,189)
tax18<-c("Staphylococcus_epidermidis",0,0,0,189)
tax19<-c("Streptococcus_agalactiae",0,0,0,189)
tax20<-c("Streptococcus_mutans",0,0,0,189)
tax21<-c("Streptococcus_pneumoniae",0,0,0,189)
otu<-rbind(otu,tax1,tax2,tax3,tax4,tax5,tax6,tax7,tax8,tax9,tax10,tax11,tax12,tax13,tax14,tax15,tax16,tax17,tax18,tax19,tax20,tax21)


#Table de taxonomie
tax<-merge(otu_seq_id[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","sequence")],id_mocks_taxsplit[,c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")],by.x = c("OTUID"),all.x=T, all.y=F)
tax<-tax[,c("sequence","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
tax$sequence<-as.character(tax$sequence)
taxo1<-c("Acinetobacter_baumannii","Bacteria","Proteobacteria","Gammaproteobacteria","Pseudomonadales","Moraxellaceae","Acinetobacter","Acinetobacter_baumannii")
taxo2<-c("Actinomyces_odontolyticus","Bacteria","Actinobacteria","Actinobacteria","Actinomycetales","Actinomycetaceae","Actinomyces","Actinomyces_odontolyticus")
taxo3<-c("Bacillus_cereus","Bacteria","Firmicutes","Bacilli","Bacillales","Bacillaceae","Bacillus","Bacillus_cereus")
taxo4<-c("Bacteroides_vulgatus","Bacteria","Bacteroidetes","Bacteroidia","Bacteroidales","Bacteroidaceae","Bacteroides","Bacteroides_vulgatus")
taxo5<-c("Clostridium_beijerinckii","Bacteria","Firmicutes","Clostridia","Clostridiales","Clostridiaceae","Clostridium","Clostridium_beijerinckii")
taxo6<-c("Deinococcus_radiodurans","Bacteria","Deinococcus-Thermus","Deinococci","Deinococcales","Deinococcaceae","Deinococcus","Deinococcus_radiodurans")
taxo7<-c("Enterococcus_faecalis","Bacteria","Firmicutes","Bacilli","Lactobacillales","Enterococcaceae","Enterococcus","Enterococcus_faecalis")
taxo8<-c("Escherichia_coli","Bacteria","Proteobacteria","Gammaproteobacteria","Enterobacteriales","Enterobacteriaceae","Escherichia","Escherichia_coli")
taxo9<-c("Helicobacter_pylori","Bacteria","Epsilonbacteraeota","Campylobacteria","Campylobacterales","Helicobacteraceae","Helicobacter","Helicobacter_pylori")
taxo10<-c("Lactobacillus_gasseri","Bacteria","Firmicutes","Bacilli","Lactobacillales","Lactobacillaceae","Lactobacillus","Lactobacillus_gasseri")
taxo11<-c("Listeria_monocytogenes","Bacteria","Firmicutes","Bacilli","Bacillales","Listeriaceae","Listeria","Listeria_monocytogenes")
taxo12<-c("Neisseria_meningitidis","Bacteria","Proteobacteria","Gammaproteobacteria","Betaproteobacteriales","Neisseriaceae","Neisseria","Neisseria_meningitidis")
taxo13<-c("Porphyromonas_gingivalis","Bacteria","Bacteroidetes","Bacteroidia","Bacteroidales","Bacteroidaceae","Bacteroides","Porphyromonas_gingivalis")
taxo14<-c("Propionibacterium_acnes","Bacteria","Actinobacteria","Actinobacteria","Propionibacteriales","Propionibacteriaceae","Cutibacterium","Propionibacterium_acnes")
taxo15<-c("Pseudomonas_aeruginosa","Bacteria","Proteobacteria","Gammaproteobacteria","Pseudomonadales","Pseudomonadaceae","Pseudomonas","Pseudomonas_aeruginosa")
taxo16<-c("Rhodobacter_sphaeroides","Bacteria","Proteobacteria","Alphaproteobacteria","Rhodobacterales","Rhodobacteraceae","Rhodobacter","Rhodobacter_sphaeroides")
taxo17<-c("Staphylococcus_aureus","Bacteria","Firmicutes","Bacilli","Bacillales","Staphylococcaceae","Staphylococcus","Staphylococcus_aureus")
taxo18<-c("Staphylococcus_epidermidis","Bacteria","Firmicutes","Bacilli","Bacillales","Staphylococcaceae","Staphylococcus","Staphylococcus_epidermidis")
taxo19<-c("Streptococcus_agalactiae","Bacteria","Firmicutes","Bacilli","Lactobacillales","Streptococcaceae","Streptococcus","Streptococcus_agalactiae")
taxo20<-c("Streptococcus_mutans","Bacteria","Firmicutes","Bacilli","Lactobacillales","Streptococcaceae","Streptococcus","Streptococcus_mutans")
taxo21<-c("Streptococcus_pneumoniae","Bacteria","Firmicutes","Bacilli","Lactobacillales","Streptococcaceae","Streptococcus","Streptococcus_pneumoniae")
tax<-rbind(tax,taxo1,taxo2,taxo3,taxo4,taxo5,taxo6,taxo7,taxo8,taxo9,taxo10,taxo11,taxo12,taxo13,taxo14,taxo15,taxo16,taxo17,taxo18,taxo19,taxo20,taxo21)

#Résumé avec les séquences en OTUID
res<-merge(otu_seq_id[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","sequence")],id_mocks_taxsplit[,c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")],by.x = c("OTUID"),all.x=T, all.y=F)
res$True_OTU<-0
res<-res[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2","True_OTU","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
res$sequence<-as.character(res$sequence)
res$Order<-as.character(res$Order)
res$Family<-as.character(res$Family)
res$Genus<-as.character(res$Genus)
res$Species<-as.character(res$Species)
res<-rbind(res,mock1,mock2,mock3,mock4,mock5,mock6,mock7,mock8,mock9,mock10,mock11,mock12,mock13,mock14,mock15,mock16,mock17,mock18,mock19,mock20,mock21)


#Export in .tsv file
#Table récapitulative
write.table(x = res, file = "/home/galati/Téléchargements/table_mocks_16S_tot.tsv",sep="\t",dec=",",row.names = F)
#Table OTU
write.table(x = otu, file = "/home/galati/Téléchargements/table_mocks_16S_otu.tsv",sep="\t",dec=",",row.names = F)
#Table Mocks
write.table(x = tax, file = "/home/galati/Téléchargements/table_mocks_16S_tax.tsv",sep="\t",dec=",",row.names = F)


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
