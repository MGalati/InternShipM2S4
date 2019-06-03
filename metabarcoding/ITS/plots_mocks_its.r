# GALATI Mathias avec l'aide de Hans Schriek
# phyloseq 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("dada2")
BiocManager::install("DESeq2")

require(phyloseq); packageVersion("phyloseq")
require(tidyverse); packageVersion("phyloseq")
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

path = "/home/galati/Téléchargements/stats/export_ITS_vsearch/export"
setwd(path)


# CSV IMPORT
otu <-  read.table(file = 'table_mocks_ITS_otu.tsv', sep = '\t', dec=".", header = TRUE, row.names = 1)
nrow(otu)
ncol(otu)
taxa <- read.table(file = 'table_mocks_ITS_tax.tsv', sep = '\t', dec=".", header = TRUE, row.names = 1)
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

plot_richness(ps)

ps.csp = subset_taxa(ps, Phylum == "Ascomycota")
plot_bar(ps.csp, fill="Genus")

gpt <- subset_taxa(ps, Kingdom=="Fungi")
gpt <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:3]), ps)
plot_heatmap(gpt, sample.label="SampleType")

plot_ordination(ps, ordinate(ps, method ="PCoS", distance = "euclidean"), color = "species") +
  geom_point(size = 3) +
  ggtitle("PCoA (log normalization)")

ordu = ordinate(ps, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps, ordu, color="sample", shape="species")


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

#BOF# Alpha div
plot_richness(ps)

#BOF# Beta div
GP.ord <- ordinate(ps, "NMDS", "bray")
p1 = plot_ordination(ps, GP.ord, type="taxa", color="Species", title="taxa")
print(p1)
p1 + facet_wrap(~Species, 3)

#NO# Plot bar / Taxonomy
ps.csp = subset_taxa(ps, Kingdom == "Bacteria")
plot_bar(ps, fill="Species")
