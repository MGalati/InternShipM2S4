# GALATI Mathias aidé par Hans Schrieke
# phyloseq 


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("dada2")
BiocManager::install("DESeq2")

require(phyloseq); packageVersion("phyloseq")
require(tidyverse); packageVersion("tidyverse")
require(reshape2); packageVersion("reshape2")
require(gridExtra); packageVersion("gridExtra")
require(scales); packageVersion("scales")
require(parallel); packageVersion("parallel")
require(permute); packageVersion("permute")
require(lattice); packageVersion("lattice")
require("plyr"); packageVersion("plyr")
require("DESeq2"); packageVersion("DESeq2")


#Note : de otu ITS à enlever et de metadata blank_8, T-16S-pomme, T-16S-pomme2, T-ITS-mangue, T-ITS-mangue2
#Vérifier à avoir biuen des _

##############################
########### DATA ############
##############################

path = "/home/galati/Téléchargements/stats/export_16S_vsearch/export"
setwd(path)


# CSV IMPORT
otu <-  read.table(file = 'table_mocks_16S_otu.tsv', sep = '\t', dec=".", header = TRUE, row.names = 1)
nrow(otu)
ncol(otu)
taxa <- read.table(file = 'table_mocks_16S_tax.tsv', sep = '\t', dec=".", header = TRUE, row.names = 1)
nrow(taxa)
ncol(taxa)

##############################
## PHYLOSEQ OBJECT CREATION ##
##############################

require(dada2); packageVersion("dada2")
require(phyloseq); packageVersion("phyloseq")
require(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw()) #set ggplot2 graphic theme 

otu <- as.matrix(otu)
otu <-otu+1
taxa <- as.matrix(taxa)
OTU = otu_table(otu, taxa_are_rows =TRUE)

row.names(taxa)
TAX = tax_table(taxa)

#Check names, they should be same
taxa_names(OTU)
taxa_names(TAX)

#Objet phyloseq
ps <- phyloseq(OTU, TAX)


##############################
######## PHYLOSEQ MANIP#######
##############################

#OK# Profondeur de séquencage
# Dessine le graphique de la distribution de séquences par OTU et par échantillon
# Crée un dataframe avec nreads : trie le nombre de reads par OTU
readsumsdf <- data.frame(nreads = sort(taxa_sums(ps), TRUE),
                         sorted = 1:ntaxa(ps), 
                         type = "OTU")

# Dessine le graphique avec ggplot
ggplot(readsumsdf, 
       aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity") + 
  scale_y_log10() 

# Crée un autre dataframe avec la profondeur de séquencage par échantillon avec sample sample_sums()
readsumsdf2 <- data.frame(nreads = sort(sample_sums(ps), TRUE), 
                          sorted = 1:nsamples(ps), 
                          type = "Samples")

# Rassemblement des deux dataframe
readsumsdf3 <- rbind(readsumsdf,readsumsdf2)

# Graph avec ggplot les données couvrent les données Type (OTU et Samples )
p  <-  ggplot(readsumsdf3, 
              aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity")
p + ggtitle("Total number of reads before Preprocessing") + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#OK# Alpha div
#Il faut appeler OTU "Formal class otu_table"
otu <-  read.table(file = 'table_mocks_16S_otu.tsv', sep = '\t', dec=".", header = TRUE, row.names = 1)
OTU = otu_table(otu, taxa_are_rows =TRUE)
estimate_richness(OTU, split = TRUE, measures = NULL)
plot_richness(OTU,measures = c("Observed", "Chao1", "ACE", "Shannon"))

##############################
######## NORMALIZATION########
##############################

# %
class(OTU)
count_tab_ps <- as.data.frame(OTU)

otu_percent_ps <- cbind(0)
for (i in 1:ncol(count_tab_ps)){
  otu_percent_ps <- cbind(otu_percent_ps,(count_tab_ps[i]/colSums(count_tab_ps[i]))*100)
  
}

sum(otu_percent_ps)
otu_percent_ps <- otu_percent_ps[,-1]

OTU = otu_table(otu_percent_ps, taxa_are_rows =TRUE)
ps <- phyloseq(OTU, TAX)

#YES# Plot bar / Taxonomy
require(ggplot2); packageVersion("ggplot2")
ps.csp = subset_taxa(ps, Kingdom == "Bacteria")
plot_bar(ps.csp, fill="Species")

p = plot_bar(ps.csp, "Sample", fill="Species")
p + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")

p = plot_bar(ps.csp, "Sample", fill="Genus")
p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

p = plot_bar(ps.csp, "Sample", fill="Phylum")
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
