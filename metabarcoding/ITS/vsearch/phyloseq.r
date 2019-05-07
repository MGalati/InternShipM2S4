# Squelette du script SCHRIEKE Hans
# Adapté par GALATI Mathias
# phyloseq 


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)

#Note : de otu ITS à enlever et de metadata blank_8, T-16S-pomme, T-16S-pomme2, T-ITS-mangue, T-ITS-mangue2
#Vérifier à avoir biuen des _

##############################
########### DATA ############
##############################

path = "/home/galati/Téléchargements/export_ITS_vsearch/export"
setwd(path)


# CSV IMPORT
otu <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE, row.names = 1)
nrow(otu)
ncol(otu)
taxa <- read.table(file = 'ITS_vsearch_taxonomy.tsv', sep = '\t', dec=".", header = TRUE, row.names = 1)
nrow(taxa)
ncol(taxa)
metadata <- read.table(file = 'metadata_ITS.tsv', sep = ',', dec=",", header = TRUE, row.names = 1)
nrow(metadata)
ncol(metadata)

# REMOVE SAMPLES WITHOUT READS

metadata <- as.data.frame(t(metadata))
otu <- otu[, order(colnames(otu))]
metadata <- metadata[, order(colnames(metadata))]

# PRINT OTU EMPTY
for (i in 1:ncol(otu)){
  if (sum(otu[,i])==0){
    print(names(otu[i]))
  }
}

#REMOVE
for (i in 1:ncol(otu)){
  if (sum(otu[,i])==0){
    print(names(otu[i]))
    otu <- otu[,-i]
    print("otu column removed")
    metadata <- metadata[,-i]
    print("metadata column removed")
  }

}

#Recap, les colonnes correspondent aux ID
ncol(otu)
ncol(metadata)

#Nouvelle transposée de metadata pour DESeq2
metadata <- t(metadata)


#metadata[metadata=="N"] <- NA # negative controls 
#metadata <- na.omit(metadata)

otu <- as.matrix(otu)
taxa <- as.matrix(taxa)
metadata <- as.data.frame(metadata)


##############################
####### NORMALIZATION ########
##############################

# Min/Max
otu_norm <- (otu -min(otu))/(max(otu)-min(otu))

# Log
otu_log <- log(otu +1)

# %
otu_percent <- (otu*100/colSums(otu))
otu_percent <- as.matrix(otu_percent)

# DESeq2
colnames(otu)
rownames(metadata)

ncol(otu)
nrow(metadata)

#otu2 <- otu[,c(1:251)] # remove the last "Undertemined" column (Pour Hans)
#otu2 <- as.matrix(otu2)

metadata2 <- metadata[order(metadata$SAMPLE),] 
metadata2 <- as.data.frame(metadata2)

rm(otu2)
otu2 <- otu
otu2 <- otu2[, order(colnames(otu2))]
otu2 <- (otu2+1) # allows to remove the zero --> DESeq doesn't work with zero 

#Vérification
ncol(otu2)
ncol(metadata2)
nrow(otu2)
nrow(metadata2)

#Check de ce qui n'est pas identique
y=0
for (i in colnames(otu)){
  x=FALSE
  for (j in rownames(metadata2)){
    if (i==j){
      x=TRUE
    }
  }
  if (x==FALSE){
    print(i)
    y=y+1
    print(y)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("XML")
install.packages("annotate")
install.packages("genefilter")
install.packages("geneplotter")
install.packages("DESeq2")
library(DESeq2)

write.table(x = otu2, file = "/home/galati/Téléchargements/otu2.csv")
any(sum(otu2)!=0)

dds <- DESeqDataSetFromMatrix(countData=otu2, colData=metadata2, design=~Exposure)
dds <- DESeq(dds)
#Erreur : impossible d'allouer un vecteur de taille 216.8 Mo

otu_deseq <- counts(dds, normalized=TRUE) 




##############################
## PHYLOSEQ OBJECT CREATION ##
##############################

library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw()) #set ggplot2 graphic theme 

otu_norm <- as.matrix(otu_norm)
otu_log <- as.matrix(otu_log)
otu_percent <- as.matrix(otu_percent)
#otu_deseq <- as.matrix(otu_deseq)

OTU = otu_table(otu, taxa_are_rows =TRUE)
OTU_norm = otu_table(otu_norm, taxa_are_rows = TRUE)
OTU_log = otu_table(otu_log, taxa_are_rows = TRUE)
OTU_percent = otu_table(otu_percent, taxa_are_rows = TRUE)
#OTU_deseq = otu_table(otu_deseq, taxa_are_rows = TRUE)

row.names(taxa)
TAX = tax_table(taxa)
SAM = sample_data(metadata)

#Check names, they should be same
taxa_names(OTU)
taxa_names(TAX)


ps <- phyloseq(OTU, TAX, SAM)
ps_norm <- phyloseq(OTU_norm, TAX, SAM) 
ps_log <- phyloseq(OTU_log, TAX, SAM)
ps_percent <- phyloseq(OTU_percent, TAX, SAM)
#ps_deseq <- phyloseq(OTU_deseq, TAX, SAM)




##############################
######### DECONTAM ###########
##############################

library(decontam); packageVersion("decontam")

summary(sample_sums(ps))
summary(taxa_sums(ps))
#summary(otu_sums(ps))

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

pdf("LibrarySize.pdf")
ggplot(data=df, aes(x=Index, y=LibrarySize, color=BLANK)) + geom_point()
dev.off()

as.numeric(get_variable(ps, "DNA"))
get_variable(ps, "DNA")
sample_data(ps)

sample_data(ps)$DNA <- as.numeric(get_variable(ps, "DNA"))


#  Identify Contaminants - FREQUENCY
contamdf.freq <- isContaminant(ps, method="frequency", conc="DNA")
head(contamdf.freq)

table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="DNA") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),3)], conc="DNA") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")


ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps)
ps.noncontam

# Identify Contaminants - PREVALENCE
any(taxa_sums(ps) == 0)
psa <- prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(psa) == 0)

sample_data(ps)$is.neg <- sample_data(ps)$BLANK == "Control sample"
contamdf.prev <- isContaminant(psa, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

contamdf.prev05 <- isContaminant(psa, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

ps.pa <- transform_sample_counts(psa, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$BLANK == "Control sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$BLANK == "True sample", ps.pa)


# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

pdf("PrevalenceDecontam.pdf")
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
dev.


# Identify Contaminants - COMBINED (frequency + prevalence)
contamdf.combi <- isContaminant(ps, method="combined", conc="DNA", neg="is.neg")
table(contamdf.combi$contaminant)

ps.noncontam <- prune_taxa(!contamdf.combi$contaminant, ps)
ps.noncontam 



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

ps.E <- subset_samples(ps.noncontam, ORGANISM == "E") # full organism 
ps.I <- subset_samples(ps.noncontam, ORGANISM == "I") # intestine
ps.GS <- subset_samples(ps.noncontam, ORGANISM == "GS") # salivary gland
ps.O <- subset_samples(ps.noncontam, ORGANISM == "O") # ovaries
ps.P <- subset_samples(ps.noncontam, ORGANISM == "P") # organs pull
ps.T <- subset_samples(ps.noncontam, LOCATION == "T") # blanks
ps.wT <- subset_samples(ps.noncontam, LOCATION != "T") # without blanks



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


##############################
####### RICHNESS PLOT ########
##############################
subset_samples(ps.noncontam, BLANK == "Control sample") %>%
  plot_richness(measures=c("Observed","Shannon","ACE")) -> toto

ggplot(toto$data)
#voir ggplot pour boxplot



# Save the session
install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("dada2.Rda")

