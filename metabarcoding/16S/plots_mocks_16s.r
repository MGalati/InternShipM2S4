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

#BOF# Alpha div
plot_richness(ps)

#BOF# Beta div
GP.ord <- ordinate(ps, "NMDS", "bray")
p1 = plot_ordination(ps, GP.ord, type="taxa", color="Phylum", title="taxa")
print(p1)
p1 + facet_wrap(~Phylum, 3)
#NO# Plot bar / Taxonomy
ps.csp = subset_taxa(ps, Kingdom == "Bacteria")
plot_bar(ps.csp, fill="Species")


# Let's explore the rarefaction curves i.e., OTU richness vs sequencing depth
ps %>%
  otu_table() %>%
  t() %>%
  vegan::rarecurve()

# We can do something nicer with ggplot
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")

ps<-as.data.frame(ps)
p <- ggrare(ps,
            step = 500,
            color = NULL,
            plot = T,
            parallel = T,
            se = F)
p <- p + 
  geom_vline(xintercept = min(sample_sums(ps)), 
             color = "gray60")
plot(p)

# La table d'OTU va être filtrée
# Exploration de la taxonomie au niveau du Kingdom
tax_table(ps)[,c("Kingdom")] %>% unique()

# Remove untargeted OTU (we consider Unclassified OTU at the Kingdom level as noise) using subset_taxa
data <- subset_taxa(ps, 
                    Kingdom != "Unclassified" &
                      Order !="Chloroplast" &
                      Family != "Mitochondria")

# Remove low occurence / abundance OTU i.e.,  more than 10 sequences in total and appearing in more than 1 sample
data <-  filter_taxa(data, 
                     function(x) sum(x >= 10) > (1), 
                     prune =  TRUE) 

# Rarefy to en even sequencing depth (i.e., min(colSums(otu_table(data)))
data_rare <- rarefy_even_depth(data, 
                               sample.size = min(colSums(otu_table(data))), 
                               rngseed = 63)

# Rarefaction curves on filtered data
p <- ggrare(data_rare, step = 50, color = "env_material", plot = T, parallel = T, se = F)
p 

# One can export the filtered OTU table
write.csv(cbind(data.frame(otu_table(data_rare)),
                tax_table(data_rare)), 
          file="filtered_otu_table.csv")

# We can now explore the alpha-dviersity on the filtered and rarefied data
p <- plot_richness(ps, 
                   x="sample", 
                   color=NULL, 
                   measures=c("Observed","Shannon","ACE"), 
                   nrow = 1)
print(p)

# That plot could be nicer
# data to plot are stored in p$data
p$data %>% head()

# boxplot using ggplot 
ggplot(p$data,aes(env_material,value,colour=env_material)) +
  facet_grid(variable ~ env_material, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=1)

# More Complex
ggplot(p$data,aes(env_material,value,colour=env_material,shape=env_material)) +
  facet_grid(variable ~ env_material, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.8, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  ylab("Diversity index")  + xlab(NULL) + theme_bw()

# Export the alpha div values into a dataframe in short format
rich.plus <- dcast(p$data,  samples + env_material ~ variable)
write.csv(rich.plus, file="alpha_div.csv")

# Alpha-div Stats using TukeyHSD on ANOVA
TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ env_material, data =  rich.plus))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$env_material)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ env_material, data =  rich.plus))))$p.value
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
                color = "env_material", 
                shape="env_material", 
                title = "PCoA sqrt Bray curtis", 
                label= "SampleID" ) + 
  geom_point(aes(size=rich.plus$Observed)) +
  theme_bw()

# Let's see if the observed pattern is significant using PERMANOVA i.e., adonis function from vegan
adonis(dist ~ get_variable(data_rare, "env_material"), permutations = 1000)$aov.tab

# Let's see if there are difference in dispersion (i.e., variance)
boxplot(betadisper(dist, 
                   get_variable(data_rare, "env_material")),las=2, 
        main=paste0("Multivariate Dispersion Test Bray-Curtis "," pvalue = ", 
                    permutest(betadisper(dist, get_variable(data_rare, "env_material")))$tab$`Pr(>F)`[1]))

# ANOSIM test can also test for differences among group 
plot(anosim(dist, get_variable(data_rare, "env_material"))
     ,main="ANOSIM Bray-Curtis "
     ,las=2)

# Now, we would like to plot the distribution of phylum transformed in %
data_rare %>% transform_sample_counts(function(x) x/sum(x)) %>%
  plot_bar(fill="Phylum") +
  facet_wrap(~ env_material, scales = "free_x", nrow = 1) +
  ggtitle("Bar plot colored by Phylum ") +
  theme(plot.title = element_text(hjust = 0.5))

# We can generate a nicer plot using plot_composition function
p <- plot_composition(data_rare,
                      taxaRank1 = "Kingdom",
                      taxaSet1 ="Bacteria",
                      taxaRank2 = "Phylum", 
                      numberOfTaxa = 20, 
                      fill= "Phylum") +
  facet_wrap(~env_material, scales = "free_x", nrow = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) 

plot(p)


##############################
########### ADONIS ###########
##############################

library(vegan)
library(dist)
dist.jac <- vegdist(otu, method="jaccard", binary=TRUE)
dist.jac <- as.matrix(dist.jac)
dist.jac<- as.data.frame(dist.jac)
dist2.jac <- distance(otu, method = "jaccard")
as.recursive(ps)
is.atomic(ps)
adonis(formula = dist.jac~LOCATION, data=metadata2, permutation = 9999)


adonis(vegdist(t(otu_table(ps)), method = "bray") ~ LOCATION,
       data=as(sample_data(ps.noncontam), "data.frame"), permutation = 9999)


# testtt <- t(otu_table(ps))
# test2 <- otu_table(ps)
# View(testtt)
# View(test2)


##############################
##### RAREFACTION CURVE ######
##############################

source('functions.R') # import amp_rarecurve, plot_composition

library(dplyr)
library(vegan)
library(magrittr)
library(ggplot2)
library(datasets)

ps.E <- subset_samples(ps, ORGANISM == "E") # full organism 
ps.I <- subset_samples(ps, ORGANISM == "I") # intestine
ps.GS <- subset_samples(ps, ORGANISM == "GS") # salivary gland
ps.O <- subset_samples(ps, ORGANISM == "O") # ovaries
ps.P <- subset_samples(ps, ORGANISM == "P") # organs pull
ps.T <- subset_samples(ps, LOCATION == "T") # blanks
ps.wT <- subset_samples(ps, LOCATION != "T") # without blanks



amp_rarecurve(subset_samples(ps, LOCATION =="N"),
              step=100,
              label = T,
              legend.position = "bottomright",
              legend = T)



readsumsdf = data.frame(nreads = sort(taxa_sums(ps.wT), TRUE), sorted = 1:ntaxa(ps.wT), 
                        type = "OTUs")

readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps.wT), TRUE), sorted = 1:nsamples(ps.wT), 
                                          type = "Samples"))

ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity") + ggtitle("Total number of reads before Preprocessing (Sans T?moins)") + scale_y_log10() +
  facet_wrap(~type, ncol = 1, scales = "free") #+ scale_y_log10()

graph2ppt(file=paste0("Preanalysis",".ppt"),append=T,width=9,aspectr=sqrt(2))

ggrare(ps.noncontam, step = 7, label = NULL, color = NULL,
       plot = TRUE, parallel = TRUE, se = TRUE)



##############################
##### COMPOSITION PLOT #######
##############################

library(scales)
library(reshape2)
library(ggplot2)

scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R",
             "NCM_fit.R")
urls <-
  paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/",
         scripts)

for (url in urls) {
  source(url)
}



subset_samples(ps.noncontam, LOCATION == "T") %>%
  plot_composition("Kingdom", "Bacteria", "Species", numberOfTaxa =
                     1000, fill = "Species")




# Save the session
install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("metab.Rda")
