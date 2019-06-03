#Préparation de données et comparaison de mocks ITS
#Mathias

####Import OTU/TAX####

# Importation ITS vsearch
path = "/home/galati/Téléchargements/stats/export_ITS_vsearch/export"
setwd(path)
otu_ITS_vsearch <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE)
taxa_ITS_vsearch <- read.table(file = 'ITS_vsearch_taxonomy.tsv', sep = '\t', dec=".", header = TRUE)

# Importation ITS deblur
path = "/home/galati/Téléchargements/stats/export_ITS_deblur/export"
setwd(path)
otu_ITS_deblur <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE)
taxa_ITS_deblur <- read.table(file = 'ITS_deblur_taxonomy.tsv', sep = '\t', dec=".", header = TRUE)

# Importation ITS dada2
path = "/home/galati/Téléchargements/stats/export_ITS_dada2/export"
setwd(path)
otu_ITS_dada2 <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE)
taxa_ITS_dada2 <- read.table(file = 'ITS_dada2_taxonomy.tsv', sep = '\t', dec=".", header = TRUE)

####Renommage####

#Renommage des mocks
colnames(otu_ITS_vsearch)[colnames(otu_ITS_vsearch)=="ITS"] <- "Mock_vsearch"
colnames(otu_ITS_deblur)[colnames(otu_ITS_deblur)=="ITS"] <- "Mock_deblur"
colnames(otu_ITS_dada2)[colnames(otu_ITS_dada2)=="ITS"] <- "Mock_dada2"

#Renommage des Taxon
names(taxa_ITS_deblur)[1]<-"OTUID"
colnames(taxa_ITS_vsearch)[colnames(taxa_ITS_vsearch)=="Taxon"] <- "Taxon_vsearch"
colnames(taxa_ITS_deblur)[colnames(taxa_ITS_deblur)=="Taxon"] <- "Taxon_deblur"
colnames(taxa_ITS_dada2)[colnames(taxa_ITS_dada2)=="Taxon"] <- "Taxon_dada2"

####Import SEQ####

library("Biostrings")

#Importation des séquences en dataframe pour les OTUID de vsearch
path = "/home/galati/Téléchargements/stats/export_ITS_vsearch/export"
setwd(path)
fastaFile <- readDNAStringSet("dna-sequences.fasta")
OTUID = names(fastaFile)
sequence_vsearch = paste(fastaFile)
s_vsearch <- data.frame(OTUID, sequence_vsearch)

#Importation des séquences en dataframe pour les OTUID de deblur
path = "/home/galati/Téléchargements/stats/export_ITS_deblur/export"
setwd(path)
fastaFile <- readDNAStringSet("dna-sequences.fasta")
OTUID = names(fastaFile)
sequence_deblur = paste(fastaFile)
s_deblur <- data.frame(OTUID, sequence_deblur)

#Importation des séquences en dataframe pour les OTUID de dada2
path = "/home/galati/Téléchargements/stats/export_ITS_dada2/export"
setwd(path)
fastaFile <- readDNAStringSet("dna-sequences.fasta")
OTUID = names(fastaFile)
sequence_dada2 = paste(fastaFile)
s_dada2 <- data.frame(OTUID, sequence_dada2)

####Merge####

#Merge des tables pour avoir comme colonnes finales OTU/Seq/Mocks
seq_merge_ITS_vsearch<-merge(s_vsearch[,c("OTUID","sequence_vsearch")],otu_ITS_vsearch[,c("OTUID","Mock_vsearch")],by.x = c("OTUID"),all.x=F, all.y=F)
seq_merge_ITS_deblur<-merge(s_deblur[,c("OTUID","sequence_deblur")],otu_ITS_deblur[,c("OTUID","Mock_deblur")],by.x = c("OTUID"),all.x=F, all.y=F)
seq_merge_ITS_dada2<-merge(s_dada2[,c("OTUID","sequence_dada2")],otu_ITS_dada2[,c("OTUID","Mock_dada2")],by.x = c("OTUID"),all.x=F, all.y=F)

seq_merge_vsearch_deblur<-merge(seq_merge_ITS_vsearch[,c("sequence_vsearch","Mock_vsearch")],seq_merge_ITS_deblur[,c("sequence_deblur","Mock_deblur")],by.x = c("sequence_vsearch"),by.y = c("sequence_deblur"),all.x=T, all.y=T)
seq_merge_all<-merge(seq_merge_vsearch_deblur[,c("sequence_vsearch","Mock_vsearch","Mock_deblur")],seq_merge_ITS_dada2[,c("sequence_dada2","Mock_dada2")],by.x = c("sequence_vsearch"),by.y = c("sequence_dada2"),all.x=T, all.y=T)

####Manip####

#Renommage de la colonne des séquences
colnames(seq_merge_all)[colnames(seq_merge_all)=="sequence_vsearch"] <- "sequence"

#Remplacement des NA avec des 0
seq_merge_all[is.na(seq_merge_all)] <- 0

#Supression des lignes où il y a que des comptes d'OTU à 0
seq_final <- seq_merge_all[!rowSums(seq_merge_all[, -1] == 0) == (ncol(seq_merge_all)-1), ]

#Merge des OTUIDs pour récupérer ensuite la taxonomie
#seq_final<-merge(seq_final[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2")],seq_merge_ITS_vsearch[,c("OTUID","sequence_vsearch")],by.x = c("sequence"),by.y = c("sequence_vsearch"),all.x=T, all.y=F)
#seq_final<-merge(seq_final[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2","OTUID")],seq_merge_ITS_deblur[,c("OTUID","sequence_deblur")],by.x = c("sequence"),by.y = c("sequence_deblur"),all.x=T, all.y=F)

####Récupération des séquences pour avoir une table combinant les séquences et les identifiants Qiime2####
colnames(seq_merge_ITS_vsearch)[colnames(seq_merge_ITS_vsearch)=="sequence_vsearch"] <- "sequence"
colnames(seq_merge_ITS_vsearch)[colnames(seq_merge_ITS_vsearch)=="Mock_vsearch"] <- "Mock"
colnames(seq_merge_ITS_deblur)[colnames(seq_merge_ITS_deblur)=="sequence_deblur"] <- "sequence"
colnames(seq_merge_ITS_deblur)[colnames(seq_merge_ITS_deblur)=="Mock_deblur"] <- "Mock"
colnames(seq_merge_ITS_dada2)[colnames(seq_merge_ITS_dada2)=="sequence_dada2"] <- "sequence"
colnames(seq_merge_ITS_dada2)[colnames(seq_merge_ITS_dada2)=="Mock_dada2"] <- "Mock"
otu_vsearch_deblur<-rbind(seq_merge_ITS_vsearch,seq_merge_ITS_deblur)
otu_vsearch_deblur_dada2<-rbind(otu_vsearch_deblur,seq_merge_ITS_dada2)

otu_seq_id<-merge(seq_final[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2")],otu_vsearch_deblur_dada2[,c("OTUID","sequence")],by.x = c("sequence"),all.x=T, all.y=F)

####Obtention d'une table OTU/Mocks/Taxonomie####
colnames(taxa_ITS_deblur)[colnames(taxa_ITS_deblur)=="Taxon_deblur"] <- "Taxon"
colnames(taxa_ITS_vsearch)[colnames(taxa_ITS_vsearch)=="Taxon_vsearch"] <- "Taxon"
colnames(taxa_ITS_vsearch)[colnames(taxa_ITS_vsearch)=="Feature.ID"] <- "OTUID"
colnames(taxa_ITS_dada2)[colnames(taxa_ITS_dada2)=="Taxon_dada2"] <- "Taxon"
colnames(taxa_ITS_dada2)[colnames(taxa_ITS_dada2)=="Feature.ID"] <- "OTUID"
tax_tot<-rbind(taxa_ITS_deblur,taxa_ITS_vsearch, taxa_ITS_dada2)
otu_mock_tax<-merge(otu_seq_id[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2")],tax_tot[,c("OTUID","Taxon")],by.x = c("OTUID"),all.x=T, all.y=T)

####Combiner les infos d'OTU et de Taxonomie####
comb_tax <- data.frame(newCol=paste(otu_mock_tax$OTUID,otu_mock_tax$Taxon,sep=";"))
names(comb_tax)[1]<-"OTUID"

####Recherche pour séparer la colonnes Taxon en différentes colonnes pour ressembler à la nommenclature du tuto####
#"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
library(stringr)
tax_split <- str_split_fixed(comb_tax$OTUID,";.__",8)
tax_split<-as.data.frame(tax_split)
names(tax_split) <- c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
#Vérification
colnames(tax_split)

id_mocks_taxsplit<-merge(tax_split[,c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")],otu_mock_tax[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2")],by.x = c("OTUID"),all.x=T, all.y=F)
id_mocks_taxsplit<-id_mocks_taxsplit[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","Kingdom","Phylum","Class","Order","Family","Genus","Species")]


####True Mock#### Tout à été en propotion égale 0.0625 cad 886651 en fonction du nombre de compte de dada2
id_mocks_taxsplit$True_OTU<-0
id_mocks_taxsplit<-id_mocks_taxsplit[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","True_OTU","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
id_mocks_taxsplit$OTUID<-as.character(id_mocks_taxsplit$OTUID)
id_mocks_taxsplit$Order<-as.character(id_mocks_taxsplit$Order)
id_mocks_taxsplit$Family<-as.character(id_mocks_taxsplit$Family)
id_mocks_taxsplit$Genus<-as.character(id_mocks_taxsplit$Genus)
id_mocks_taxsplit$Species<-as.character(id_mocks_taxsplit$Species)

mock1<-c("Candida_catenulata",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Saccharomycetales_fam_Incertae_sedis","Candida","Candida_catenulata")
mock2<-c("Debaryomyces_hansenii",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Debaryomycetaceae","Debaryomyces","Debaryomyces_hansenii")
mock3<-c("Dipodascus_geotrichum_Type_1",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Dipodascacea","Geotrichum","Dipodascus_geotrichum_Type_1")
mock4<-c("Dipodascus_geotrichum_Type_2",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Dipodascacea","Geotrichum","Dipodascus_geotrichum_Type_2")
mock5<-c("Dipodascus_geotrichum_Type_3",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Dipodascacea","Geotrichum","Dipodascus_geotrichum_Type_3")
mock6<-c("Dipodascus_geotrichum_Type_4",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Dipodascacea","Geotrichum","Dipodascus_geotrichum_Type_4")
mock7<-c("Kluyveromyces_lactis",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Saccharomycetaceae","Kluyveromyces","Kluyveromyces_lactis")
mock8<-c("Cyberlindnera_jadinii",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Phaffomycetaceae","Cyberlindnera","Cyberlindnera_jadinii")
mock9<-c("Zygosaccharomyces_rouxii",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Saccharomycetaceae","Zygosaccharomyces","Zygosaccharomyces_rouxii")
mock10<-c("Hyphopichia_burtonii",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Debaryomycetaceae","Hyphopichia","Hyphopichia_burtonii")
mock11<-c("Pichia_kudriavzevii",0,0,0,886651,"Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Pichiaceae","Pichia","Pichia_kudriavzevii")
mock12<-c("Penicillium_roqueforti",0,0,0,886651,"Fungi","Ascomycota","Eurotiomycetes","Eurotiales","Aspergillaceae","Penicillium","Penicillium_roqueforti")
mock13<-c("Penicillium_allii",0,0,0,886651,"Fungi","Ascomycota","Eurotiomycetes","Eurotiales","Aspergillaceae","Penicillium","Penicillium_allii")
mock14<-c("Penicillium_commune",0,0,0,886651,"Fungi","Ascomycota","Eurotiomycetes","Eurotiales","Aspergillaceae","Penicillium","Penicillium_commune")
mock15<-c("Fusarium_domesticum",0,0,0,886651,"Fungi","Ascomycota","Sordariomycetes","Hypocreales","Nectriaceae","Fusarium","Fusarium_domesticum")
mock16<-c("Scopulariopsis_fusca",0,0,0,886651,"Fungi","Ascomycota","Sordariomycetes","Microascales","Microascaceae","Scopulariopsis","Scopulariopsis_fusca")

id_mocks_taxsplit<-rbind(id_mocks_taxsplit,mock1,mock2,mock3,mock4,mock5,mock6,mock7,mock8,mock9,mock10,mock11,mock12,mock13,mock14,mock15,mock16)

####Création des tables pour phyloseq####
#Table OTU
otu<-otu_seq_id[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2")]
otu$sequence<-as.character(otu$sequence)
otu$True_OTU<-0
tax1<-c("Candida_catenulata",0,0,0,886651)
tax2<-c("Debaryomyces_hansenii",0,0,0,886651)
tax3<-c("Dipodascus_geotrichum_Type_1",0,0,0,886651)
tax4<-c("Dipodascus_geotrichum_Type_2",0,0,0,886651)
tax5<-c("Dipodascus_geotrichum_Type_3",0,0,0,886651)
tax6<-c("Dipodascus_geotrichum_Type_4",0,0,0,886651)
tax7<-c("Kluyveromyces_lactis",0,0,0,886651)
tax8<-c("Cyberlindnera_jadinii",0,0,0,886651)
tax9<-c("Zygosaccharomyces_rouxii",0,0,0,886651)
tax10<-c("Hyphopichia_burtonii",0,0,0,886651)
tax11<-c("Pichia_kudriavzevii",0,0,0,886651)
tax12<-c("Penicillium_roqueforti",0,0,0,886651)
tax13<-c("Penicillium_allii",0,0,0,886651)
tax14<-c("Penicillium_commune",0,0,0,886651)
tax15<-c("Fusarium_domesticum",0,0,0,886651)
tax16<-c("Scopulariopsis_fusca",0,0,0,886651)
otu<-rbind(otu,tax1,tax2,tax3,tax4,tax5,tax6,tax7,tax8,tax9,tax10,tax11,tax12,tax13,tax14,tax15,tax16)


#Table de taxonomie
tax<-merge(otu_seq_id[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","sequence")],id_mocks_taxsplit[,c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")],by.x = c("OTUID"),all.x=T, all.y=F)
tax<-tax[,c("sequence","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
tax$sequence<-as.character(tax$sequence)
taxo1<-c("Candida_catenulata","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Saccharomycetales_fam_Incertae_sedis","Candida","Candida_catenulata")
taxo2<-c("Debaryomyces_hansenii","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Debaryomycetaceae","Debaryomyces","Debaryomyces_hansenii")
taxo3<-c("Dipodascus_geotrichum_Type_1","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Dipodascacea","Geotrichum","Dipodascus_geotrichum_Type_1")
taxo4<-c("Dipodascus_geotrichum_Type_2","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Dipodascacea","Geotrichum","Dipodascus_geotrichum_Type_2")
taxo5<-c("Dipodascus_geotrichum_Type_3","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Dipodascacea","Geotrichum","Dipodascus_geotrichum_Type_3")
taxo6<-c("Dipodascus_geotrichum_Type_4","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Dipodascacea","Geotrichum","Dipodascus_geotrichum_Type_4")
taxo7<-c("Kluyveromyces_lactis","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Saccharomycetaceae","Kluyveromyces","Kluyveromyces_lactis")
taxo8<-c("Cyberlindnera_jadinii","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Phaffomycetaceae","Cyberlindnera","Cyberlindnera_jadinii")
taxo9<-c("Zygosaccharomyces_rouxii","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Saccharomycetaceae","Zygosaccharomyces","Zygosaccharomyces_rouxii")
taxo10<-c("Hyphopichia_burtonii","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Debaryomycetaceae","Hyphopichia","Hyphopichia_burtonii")
taxo11<-c("Pichia_kudriavzevii","Fungi","Ascomycota","Saccharomycetes","Saccharomycetales","Pichiaceae","Pichia","Pichia_kudriavzevii")
taxo12<-c("Penicillium_roqueforti","Fungi","Ascomycota","Eurotiomycetes","Eurotiales","Aspergillaceae","Penicillium","Penicillium_roqueforti")
taxo13<-c("Penicillium_allii","Fungi","Ascomycota","Eurotiomycetes","Eurotiales","Aspergillaceae","Penicillium","Penicillium_allii")
taxo14<-c("Penicillium_commune","Fungi","Ascomycota","Eurotiomycetes","Eurotiales","Aspergillaceae","Penicillium","Penicillium_commune")
taxo15<-c("Fusarium_domesticum","Fungi","Ascomycota","Sordariomycetes","Hypocreales","Nectriaceae","Fusarium","Fusarium_domesticum")
taxo16<-c("Scopulariopsis_fusca","Fungi","Ascomycota","Sordariomycetes","Microascales","Microascaceae","Scopulariopsis","Scopulariopsis_fusca")
tax<-rbind(tax,taxo1,taxo2,taxo3,taxo4,taxo5,taxo6,taxo7,taxo8,taxo9,taxo10,taxo11,taxo12,taxo13,taxo14,taxo15,taxo16)

#Résumé avec les séquences en OTUID
res<-merge(otu_seq_id[,c("OTUID","Mock_vsearch","Mock_deblur","Mock_dada2","sequence")],id_mocks_taxsplit[,c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")],by.x = c("OTUID"),all.x=T, all.y=F)
res$True_OTU<-0
res<-res[,c("sequence","Mock_vsearch","Mock_deblur","Mock_dada2","True_OTU","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
res$sequence<-as.character(res$sequence)
res$Order<-as.character(res$Order)
res$Family<-as.character(res$Family)
res$Genus<-as.character(res$Genus)
res$Species<-as.character(res$Species)
res<-rbind(res,mock1,mock2,mock3,mock4,mock5,mock6,mock7,mock8,mock9,mock10,mock11,mock12,mock13,mock14,mock15,mock16)


#Export in .tsv file
#Table récapitulative
res<-as.data.frame(res)
write.table(x = res, file = "/home/galati/Téléchargements/stats/export_ITS_vsearch/export/table_mocks_ITS_tot.tsv",sep="\t",dec=",",row.names = F)
#Table OTU
otu<-as.data.frame(otu)
write.table(x = otu, file = "/home/galati/Téléchargements/stats/export_ITS_vsearch/export/table_mocks_ITS_otu.tsv",sep="\t",dec=",",row.names = F)
#Table Mocks
tax<-as.data.frame(tax)
write.table(x = tax, file = "/home/galati/Téléchargements/stats/export_ITS_vsearch/export/table_mocks_ITS_tax.tsv",sep="\t",dec=",",row.names = F)


#Visualisation des données attendues
library(phyloseq)
data(GlobalPatterns)
colnames(tax_table(GlobalPatterns))


"______________________________________________________________________________________________________________________________________________________"
