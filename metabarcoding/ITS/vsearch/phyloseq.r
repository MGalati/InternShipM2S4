# Squelette du script SCHRIEKE Hans
# Adapté par GALATI Mathias
# phyloseq 

source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)

##############################
####### NORMALIZATION ########
##############################

path = "/home/galati/Téléchargements/export_ITS_vsearch/export"
setwd(path)


# CSV IMPORT
otu <-  read.table(file = 'ASV-table.biom.tsv', sep = '\t', dec=".", header = TRUE, row.names = 1)
taxa <- read.table(file = 'ITS_vsearch_taxonomy.tsv', sep = '\t', dec=".", header = TRUE)
metadata <- read.table(file = 'metadata.tsv', sep = '\t', dec=",", header = TRUE, row.names = 1)


# REMOVE SAMPLES WITHOUT READS
metadata <- as.data.frame(t(metadata))
otu <- otu[, order(colnames(otu))]
metadata <- metadata[, order(colnames(metadata))]

# PRINT OTU EMPTY
#otu[, 1] <- as.numeric(as.character(otu[, 1]))
for (i in 1:ncol(otu)){
  if (sum(otu[,i])==0){
    print(names(otu[i]))
  }
}

#REMOVE EMPTY OTU
for (i in 1:ncol(otu)){
  if (sum(otu[,i])==0){
    print(names(otu[i]))
    otu <- otu[,-i]
    print("otu column removed")
    metadata <- metadata[,-i]
    print("metadata column removed")
  }
}

metadata <- t(metadata)

ncol(otu)
nrow(metadata)

#metadata[metadata=="N"] <- NA # negative controls 
#metadata <- na.omit(metadata)

otu <- as.matrix(otu)
taxa <- as.matrix(taxa)
metadata <- as.data.frame(metadata)

# NORMALIZATION MIN MAX
otu_norm <- (otu -min(otu))/(max(otu)-min(otu))

# NORMALIZATION LOG
otu_log <- log(otu +1)

# NORMALIZATION % 
otu_percent <- (otu*100/colSums(otu))
otu_percent <- as.matrix(otu_percent)


# NORMALIZATION DESEQ2
rownames(metadata)
colnames(otu)

otu2 <- otu[,c(1:27)] # remove "Undertemined" column
otu2 <- as.matrix(otu2)

metadata2 <- metadata[order(metadata$SAMPLE),]
metadata2 <- as.data.frame(metadata2)

otu2 <- otu2[, order(colnames(otu2))]
otu2 <- (otu2+1) #allow to remove the zero --> DESeq doesn't work with zero 

nrow(otu2)
ncol(otu2)
nrow(metadata2)

#Inversion row/col
otu3 <- t(otu2)
ncol(otu3)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
ncol(otu)
nrow(metadata2)
metadata3 <- t(metadata2)
ncol(metadata3)
dds <- DESeqDataSetFromMatrix(countData=otu, colData=metadata3, design=~Exposure)
dds <- DESeq(dds)
otu_deseq <- counts(dds, normalized=TRUE)


# Phyloseq object creation
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw()) #set ggplot2 graphic theme 

otu_norm <- as.matrix(otu_norm)
otu_log <- as.matrix(otu_log)
otu_percent <- as.matrix(otu_percent)
otu_deseq <- as.matrix(otu_deseq)

OTU = otu_table(otu, taxa_are_rows =TRUE)
OTU_norm = otu_table(otu_norm, taxa_are_rows = TRUE)
OTU_log = otu_table(otu_log, taxa_are_rows = TRUE)
OTU_percent = otu_table(otu_percent, taxa_are_rows = TRUE)
OTU_deseq = otu_table(otu_deseq, taxa_are_rows = TRUE)

TAX = tax_table(taxa)
SAM = sample_data(metadata)

ps <- phyloseq(OTU, TAX, SAM)
ps_norm <- phyloseq(OTU_norm, TAX, SAM) 
ps_log <- phyloseq(OTU_log, TAX, SAM)
ps_percent <- phyloseq(OTU_percent, TAX, SAM)
ps_deseq <- phyloseq(OTU_deseq, TAX, SAM)



# Adonis -- doesn't work for now
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
       data=as(sample_data(ps), "data.frame"), permutation = 9999)

testtt <- t(otu_table(ps))
test2 <- otu_table(ps)
View(testtt)
View(test2)



# RAREFACTION CURVE

source('functions.R') # import amp_rarecurve, plot_composition

library(dplyr)
library(vegan)
library(magrittr)

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

ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity") + ggtitle("Total number of reads before Preprocessing (Sans Témoins)") + scale_y_log10() +
  facet_wrap(~type, ncol = 1, scales = "free") #+ scale_y_log10()

graph2ppt(file=paste0("Preanalysis",".ppt"),append=T,width=9,aspectr=sqrt(2))


# DECONTAM 
library(decontam); packageVersion("decontam")

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=BLANK)) + geom_point()

summary(sample_sums(ps))
summary(taxa_sums(ps))

as.numeric(get_variable(ps, "DNA"))
get_variable(ps, "DNA")
sample_data(ps)

sample_data(ps)$DNA <- as.numeric(get_variable(ps, "DNA"))


#  Identify Contaminants - Frequency
contamdf.freq <- isContaminant(ps, method="frequency", conc="DNA")
head(contamdf.freq)

table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="DNA") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),3)], conc="DNA") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

sample_data(ps)$is.neg <- sample_data(ps)$BLANK == "Control sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps)
ps.noncontam

# Make phyloseq object of presence-absence in negative controls and true samples

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

ps.pa.neg <- prune_samples(sample_data(ps.pa)$BLANK == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$BLANK == "True Sample", ps.pa)



# PLOT COMPOSITION 

library(scales)
library(reshape2)
library(ggplot2)

subset_samples(ps, LOCATION == "T") %>%
  plot_composition("Kingdom", "Bacteria", "Species", numberOfTaxa =
                     1000, fill = "Species")

subset_samples(ps, BLANK == "Control sample") %>%
  plot_richness(measures=c("Observed","Shannon","ACE")) -> toto

ggplot(toto$data)
#voir ggplot pour boxplot



ggrare(ps, step = 7, label = NULL, color = NULL,
       plot = TRUE, parallel = FALSE, se = TRUE)

# Save the session
install.packages("session", repos = "http://cran.us.r-project.org")
library(session); packageVersion("session")
save.session("dada2.Rda")
