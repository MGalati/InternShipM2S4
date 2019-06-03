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

# Importation ITS vsearch
path = "/home/galati/Téléchargements/stats/export_16S_dada2/export"
setwd(path)
otu <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE, row.names=1)
tax <- read.table(file = '16S_dada2_taxonomy.tsv', sep = '\t', dec=".", header = TRUE)
meta <- read.table(file = 'metadata_16S.tsv', sep = '\t', dec=".", header = TRUE, row.names=1)
tree <- read_tree("rooted-tree.nwk")

#A lancer 2 fois car il y a 2 points
colnames(otu) <- sub("\\.","_", colnames(otu))

##############################
### TRAITEMENT TABLE TAXO ###
##############################
####Combiner les infos d'OTU et de Taxonomie####
tax <- data.frame(newCol=paste(tax$OTUID,tax$Taxon,sep=";"))
names(tax)[1]<-"OTUID"

####Recherche pour séparer la colonnes Taxon en différentes colonnes pour ressembler à la nommenclature du tuto####
#"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
library(stringr)
taxon <- str_split_fixed(tax$OTUID,";*D_.__",8)
taxon<-as.data.frame(taxon)
names(taxon) <- c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")

##############################
## PHYLOSEQ OBJECT CREATION ##
##############################

require(dada2); packageVersion("dada2")
require(phyloseq); packageVersion("phyloseq")
require(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw()) #set ggplot2 graphic theme 

otu <- as.matrix(otu)
otu <- (otu+1)
taxon <- as.matrix(taxon)
row.names(taxon) <- taxon[,1]
names(taxon) <- c("OTUID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxon<-taxon[,c("Kingdom","Phylum","Class","Order","Family","Genus","Species")]


OTU = otu_table(otu, taxa_are_rows =TRUE)

row.names(taxon)
TAX = tax_table(taxon)
SAM = sample_data(meta)

#Check names, they should be same
taxa_names(OTU)
taxa_names(TAX)

sample_names(OTU)
sample_names(SAM)

#Objet phyloseq
ps <- phyloseq(OTU, TAX, SAM, tree)


plot_richness(ps, x="Tree_type", measures=c("Observed","Chao1","Shannon"), color="Model")

ps.csp = subset_taxa(ps, Kingdom == "Bacteria")
plot_bar(ps, "Tree_type", fill="Phylum",facet_grid=~Sample)

gpt <- subset_taxa(ps, Kingdom=="Bacteria")
gpt <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:3]), ps)
plot_heatmap(ps, sample.label="SampleType")

plot_ordination(ps, ordinate(ps, method ="PCoS", distance = "euclidean"), color = "species") +
  geom_point(size = 3) +
  ggtitle("PCoA (log normalization)")

ordu = ordinate(ps, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps, ordu, color="Model", shape="Tree_type")

p4 = plot_ordination(ps, ordu, type="split", color="Model", shape="Tree_type", label="SampleType", title="split") 
p4


##############################
######## PHYLOSEQ MANIP#######
##############################
# ps est un objet phyloseq 
ps

# Affiche les premières lignes de ps
head(otu_table(ps))

# Il est possible de formater autrement avec le package tidy
ps %>% 
  otu_table() %>%
  head()

# Affiche la profondeur de séquencage de ps
ps %>% 
  otu_table() %>%
  colSums()

# Somme les comptages d'OTUs de chaque colonne
ps %>% sample_sums


# Crée un graphique de la profondeur de séquencage 
ps %>% 
  otu_table() %>%
  colSums() %>%
  sort() %>% 
  barplot(las=2)

# Combien de reads représentent les 10 premiers OTUs
sort(rowSums(otu_table(ps)), decreasing = T)[1:10]

# Affiche la taxonomie des premières lignes
ps %>% 
  tax_table() %>%
  head()


# Les métadonnées sont aussi stockées dans ps
sample_data(ps)$Tree_type


# Phyloseq possède quelques focntions
rank_names(ps) # Niveaux taxonomique
nsamples(ps) # Nombre d'échantillons
ntaxa(ps) # Nombre d'OTUs
sample_variables(ps) # Métadonnées

# Dessine le graphique de la distribution de séquences par OTU et par échantillon
# Crée un dataframe avec nreads : trie le nombre de reads par OTU
readsumsdf <- data.frame(nreads = sort(taxa_sums(ps), TRUE),
                         sorted = 1:ntaxa(ps), 
                         type = "OTU")

# Première ligne du dataframe
readsumsdf %>% head()


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

# Premières lignes
readsumsdf3 %>% head()

# Dernières lignes
readsumsdf3 %>% tail()

# Graph avec ggplot les données couvrent les données Type (OTU et Samples )
p  <-  ggplot(readsumsdf3, 
              aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity")
p + ggtitle("Total number of reads before Preprocessing") + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

# Let's explore the rarefaction curves i.e., OTU richness vs sequencing depth
#ps %>%
#  otu_table() %>%
#  t() %>%
#  vegan::rarecurve()

ps %>% 
  tax_table() %>%
  head()


# We can do something nicer with ggplot
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")

p <- ggrare(ps,
            step = 10,
            color = "Tree_type",
            plot = T,
            parallel = T,
            se = F)
p <- p + 
  facet_wrap(~ Tree_type ) + 
  geom_vline(xintercept = min(sample_sums(ps)), 
             color = "gray60")
plot(p)

# We are now going to filter the OTU table
# Explore the Taxonomy at the Kingdom level
tax_table(ps)[,c("Kingdom")] %>% unique()


# Remove untargeted OTU (we consider Unclassified OTU at the Kingdom level as noise) using subset_taxa
data <- subset_taxa(ps, 
                    Kingdom != "Unclassified" &
                      Order !="Chloroplast" &
                      Family != "Mitochondria")

# Remove low occurence / abundance OTU i.e.,  more than 10 sequences in total and appearing in more than 1 sample
data <-  filter_taxa(ps, 
                     function(x) sum(x >= 10) > (1), 
                     prune =  TRUE) 

# Rarefy to en even sequencing depth (i.e., min(colSums(otu_table(data)))
data_rare <- rarefy_even_depth(ps, 
                               sample.size = min(colSums(otu_table(ps))), 
                               rngseed = 63)

# Rarefaction curves on filtered data
p <- ggrare(data_rare, step = 50, color = "Tree_type", plot = T, parallel = T, se = F)
p 

# One can export the filtered OTU table
write.csv(cbind(data.frame(otu_table(data_rare)),
                tax_table(data_rare)), 
          file="filtered_otu_table.csv")


# We can now explore the alpha-dviersity on the filtered and rarefied data
p <- plot_richness(data_rare, 
                   x="sample", 
                   color="Tree_type", 
                   measures=c("Observed","Shannon","ACE"), 
                   nrow = 1)
print(p)


# That plot could be nicer

# data to plot are stored in p$data
p$data %>% head()


# boxplot using ggplot 
ggplot(p$data,aes(Tree_type,value,colour=Tree_type)) +
  facet_grid(variable ~ Tree_type, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=1)

# More Complex
ggplot(p$data,aes(Tree_type,value,colour=Tree_type,shape=Tree_type)) +
  facet_grid(variable ~ Tree_type, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  ylab("Diversity index")  + xlab(NULL) + theme_bw()

# Export the alpha div values into a dataframe in short format
rich.plus <- dcast(p$data,  samples + Tree_type ~ variable)
write.csv(rich.plus, file="alpha_div.csv")


# Alpha-div Stats using TukeyHSD on ANOVA
TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ Tree_type, data =  rich.plus))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$Tree_type)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ Tree_type, data =  rich.plus))))$p.value
TukeyHSD_Observed_df


# beta-diversity
# Compute dissimilarity
data_rare %>% transform_sample_counts(function(x) x/sum(x)) %>%
  otu_table() %>%
  t() %>%
  sqrt() %>%
  as.data.frame() %>%
  vegdist(binary=F, method = "bray") -> dist


# run PCoA ordination on the generated distance
ord <- ordinate(data_rare,"PCoA",dist)


# samples coordinate on the PCoA vecotrs are stored in but plot_ordination can make use of ord object easily
ord$vectors

plot_ordination(data_rare, 
                ord,
                color = "Tree_type", 
                shape="Tree_type", 
                title = "PCoA sqrt Bray curtis", 
                label= "SampleID" ) + 
  geom_point(aes(size=rich.plus$Observed)) +
  theme_bw()


# Let's see if the observed pattern is significant using PERMANOVA i.e., adonis function from vegan
adonis(dist ~ get_variable(data_rare, "Tree_type"), permutations = 1000)$aov.tab

# Let's see if there are difference in dispersion (i.e., variance)
boxplot(betadisper(dist, 
                   get_variable(data_rare, "Tree_type")),las=2, 
        main=paste0("Multivariate Dispersion Test Bray-Curtis "," pvalue = ", 
                    permutest(betadisper(dist, get_variable(data_rare, "Tree_type")))$tab$`Pr(>F)`[1]))

# ANOSIM test can also test for differences among group 
plot(anosim(dist, get_variable(data_rare, "Tree_type"))
     ,main="ANOSIM Bray-Curtis "
     ,las=2)

# Now, we would like to plot the distribution of phylum transformed in %
data_rare %>% transform_sample_counts(function(x) x/sum(x)) %>%
  plot_bar(fill="Phylum") +
  facet_wrap(~ Tree_type, scales = "free_x", nrow = 4) +
  ggtitle("Bar plot colored by Phylum ") +
  theme(plot.title = element_text(hjust = 0.5))

# We can generate a nicer plot using plot_composition function
p <- plot_composition(data_rare,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~Tree_type, scales = "free_x", nrow = 4) + 
  theme(plot.title = element_text(hjust = 0.5)) 

plot(p)
