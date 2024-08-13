                ### alpha diversity for otu table singleton edited

library(vegan)
library(phyloseq)
library(tidyverse)
library(patchwork)
library(agricolae)
library(rcompanion)
library(ggplot2)
library(writexl)

#read otu, metadata and taxonomy
data_grp <- read.table("C:/Users/Princesa/Desktop/meta2 (1).txt", header=TRUE, stringsAsFactors = TRUE,row.names = 1)
#row.names = 1 so that the first row is not considerate data but header

data_otu <- read.table("C:/Users/Princesa/Desktop/output_sing_editedfinal (1).txt", header = TRUE,row.names =1 )

data_otu1 <- as.data.frame(t(data_otu)) 


data_taxo <- read.table("C:/Users/Princesa/Desktop/sing_edited_tax_final (1).txt", header = TRUE, row.names = 1)


#create the phylosec objets with our three tables

OTU = otu_table(as.matrix(data_otu1), taxa_are_rows = FALSE)
# create the occurrence table object in phyloseq format

SAM = sample_data(data_grp, errorIfNULL = TRUE)
# create the sample metadata object in phyloseq format

TAX = tax_table(as.matrix(data_taxo))
## create the observation metadata object (OTU taxonomy) in phyloseq format


data_phylo <- phyloseq(OTU, SAM, TAX)
# create the phyloseq object including occurrence table data and sample/observation metadata

data_phylo
# print information about the phyloseq object

##########################################################################################
#Global exploration


library(ggfortify)
mtcars.pca <- prcomp(data_grp[,c(8,9)], center = TRUE,scale. = TRUE)
summary(mtcars.pca)

autoplot(mtcars.pca, data=data_grp, 
         label.size = 3,
         colour='siteic',
         label= TRUE, loadings = TRUE,
         loadings.label = TRUE)

sorted <- data_grp[order(data_grp$siteic, decreasing= TRUE),]

mtcars.pca <- prcomp(sorted[c(1:16),c(8,9)], center = TRUE,scale. = TRUE)
summary(mtcars.pca)

autoplot(mtcars.pca, data=sorted[c(1:16),], 
         label.size = 3,
         colour='siteic',
         label= TRUE, loadings = TRUE,
         loadings.label = TRUE)

mtcars.pca <- prcomp(sorted[c(17:40),c(8,9)], center = TRUE,scale. = TRUE)
summary(mtcars.pca)

autoplot(mtcars.pca, data=sorted[c(17:40),], 
         label.size = 3,
         colour='siteic',
         label= TRUE, loadings = TRUE,
         loadings.label = TRUE)

#OTU table

sum(data_otu1) 
#total reads on the OTU table = 1862697

#How many samples and variables there are in the OTU table:
nb_samples <- dim(data_otu1)[1] # nb of columns, here samples 
nb_samples
#Results: 40

nb_var <- dim(data_otu1)[2] # number of rows, here variables (OTU)
nb_var
#Results: 1578

#Sample metadata table

#samples and factors are in the metadata table:
dim(data_grp)[1] # number of samples. 
nb_factors <- dim(data_grp)[2] # number of factors
nb_factors
#Samples: 40
#Factors: 4

#How many samples per treatment:
summary(data_grp)



#OTU data properties

#Number of zeros and percentage of zeros in the OTU table
sum(data_otu1 == 0)
sum(data_otu1 == 0) / (nb_var * nb_samples) * 100

#A zero value in the OTU table means that we could not sequence any read for a specific OTU in a given sample


#Visualize the count frequency in the OTU table using a histogram.

hist(as.matrix(data_otu1), 
     max(data_otu1), 
     right = FALSE, 
     las = 1, 
     xlab = "Occurrence value", 
     ylab = "Frequency", 
     main = "Occurrence frequency",
     col="darkmagenta")

h <- hist(as.matrix(data_otu1), 
          max(data_otu1), 
          right = FALSE, 
          las = 1, 
          xlab = "Occurrence value", 
          ylab = "Frequency", 
          main = "Occurrence frequency",
          col="darkmagenta")
h


#Number of non zero values per OTU

#In order to check how the different OTU are shared bewteen samples_

non_zero <- 0*1:nb_var

for (i in 1:nb_var){
  non_zero[i]<-sum(data_otu1[,i] != 0)
}

plot(non_zero, xlab = "OTU", ylab = "Frequency", main = "Number of non zero values", las = 1)


############################################################################################
#Sequencing depth

#Sequencing depth and library size represent the number of counts per sample
#The rarefaction curve shows how many new OTU are observed when we obtain new reads
#for a given sample. If the sequencing depth is enough, we should observe a plateau,
#meaning that even if we sequence new reads they will belong to OTUs already observed
rarecurve((data_otu1[c(1,2,3,4),]), step = 1000, cex = 0.75, las = 1,main="PAE sintom?tico")
rarecurve((data_otu1[c(5,6,7,8),]), step = 1000, cex = 0.75, las = 1,main="TG controlo")
rarecurve((data_otu1[c(9,10,11,12),]), step = 1000, cex = 0.75, las = 1,main="PG sintom?tico")
rarecurve((data_otu1[c(13,14,15,16),]), step = 1000, cex = 0.75, las = 1,main="PGD controlo")
rarecurve((data_otu1[c(17,18,19,20),]), step = 1000, cex = 0.75, las = 1,main="PG controlo")
rarecurve((data_otu1[c(21,22,23,24),]), step = 1000, cex = 0.75, las = 1,main="TG sintom?tico")
rarecurve((data_otu1[c(25,26,27,28),]), step = 1000, cex = 0.75, las = 1,main="Cool controlo")
rarecurve((data_otu1[c(29,30,31,32),]), step = 1000, cex = 0.75, las = 1,main="PGD sintom?tico")
rarecurve((data_otu1[c(33,34,35,36),]), step = 1000, cex = 0.75, las = 1,main="PAE controlo")



rarecurve((data_otu1[c(1,2,3,4,9,10,11,12,21 ,22,23,24,29,30,31,32,37,38,39,40),]), step = 1000, cex = 0.75, las = 1,main="Sintomáticos")
rarecurve((data_otu1[c(5,6,7,8,13,14,15,16,17,18,19,20,25,26,27,28,33,34,35,36),]), step = 1000, cex = 0.75, las = 1,main="Controlo")

rarecurve((data_otu1[c(1,2,3,4,9,10,11,12,21 ,22,23,24,29,30,31,32,37,38,39,40,5,6,7,8,13,14,15,16,17,18,19,20,25,26,27,28,33,34,35,36),]), step = 1000, cex = 0.75, las = 1,main="Sintomáticos", col="blue")

rarecurve((data_otu1[c(1,2,3,4,5,6,7,8,13,14,15,16,21,22,24,23,29,30,31,32,33,34,35,36),]), step = 1000, cex = 0.75, las = 1,main="Bafáta")
rarecurve((data_otu1[c(9,10,11,12,17,18,19,20,25,26,27,28,37,38,39),]), step = 1000, cex = 0.75, las = 1,main="Biombo")



rarecurve(data_otu1, step = 1000, cex = 0.75, las = 1)
#We will now plot the total number of counts per sample.

#sum_seq <- rowSums(data_otu1[ ,5])
#plot(sum_seq, ylim=c(0,90000), main=c("Number of counts per sample"), xlab=c("Samples"))
#sum_seq


par(cex = 1, las = 1)
lty_vector <- c(2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1)
rarecurve(data_otu1, 
          step = 100, 
          lwd=1.3, 
          xlab = "Number of reads", 
          ylab = "Number of OTU observed", 
          xlim = c(-50, 80000), 
          ylim = c(-5, 1500), 
          label = FALSE, 
          lty = lty_vector) 
#Some samples reach a plateau, so the sequencing depth was enough

#Library size

#We will now plot the total number of counts per sample:

sum_seq <- rowSums(data_otu1)
plot(sum_seq, ylim=c(0,90000), main=c("Number of counts per sample"), xlab=c("Samples"))
sum_seq

#Qual a sample com o valor máximo e mínimo
min(sum_seq)
#Sample 12 =18719
max(sum_seq)
#Sample 13 = 78253

#We reach the conclusion we need normalization



####################################################################################
#Alpha-diversity

#Alpha-diversity represents diversity within an ecosystem or a sample,
#in other words, what is there and how much is there in term of species. 
#However, it is not easy to define a species and we can calculate alpha-diversity at different taxonomic levels.
#we are looking at the OTU level (clustered at 97% similarity thresholds).

#Several alpha-diversity indices can be calculated. Within the most commonly used:

#Richness represents the number of species observed in each sample.
#Chao1 estimates the total richness.
#Pielou’s evenness provides information about the equity in species abundance in each sample,
#in other words are some species dominating others or do all species have quite the same abundances.
#Shannon index provides information about both richness and evenness.


#It is important to not use filtered data because many richness estimates are modeled on singletons and doubletons in the occurrence table.
#So, you need to leave them in the dataset if you want a meaningful estimate.


#Indices calculation

data_richness <- estimateR(data_otu1) # calculate richness and Chao1 using vegan package

data_evenness <- diversity(data_otu1) / log(specnumber(data_otu1))# calculate evenness index using vegan package

data_shannon <- diversity(data_otu1, index = "shannon") # calculate Shannon index using vegan package

data_alphadiv <- cbind(data_grp, t(data_richness), data_shannon, data_evenness) # combine all indices in one data table

rm(data_richness, data_evenness, data_shannon) 

#Putting it in a tidy format
data_alphadiv_tidy <- 
  data_alphadiv %>%
  mutate(sample_id = rownames(data_alphadiv)) %>%
  gather(key   = alphadiv_index,
         value = obs_values,
         -ID, -Date, -Tabanca, -Site, -ic)

head(data_alphadiv_tidy)
#This give us all the values for the indexes we calculated by sample

#Visualization, graphics
#by site

P1 <- ggplot(data_alphadiv, aes(x=ic, y=S.obs)) +
  geom_boxplot(fill=c("blue","red")) +
  labs(title= 'Richness', x= ' ', y= '', tag = "A") +
  geom_point()


P2 <- ggplot(data_alphadiv, aes(x=ic, y=S.chao1)) +
  geom_boxplot(fill=c("blue","red")) +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  geom_point()

P3 <- ggplot(data_alphadiv, aes(x=ic, y=data_evenness)) +
  geom_boxplot(fill=c("blue","red")) +
  labs(title= 'Eveness', x= ' ', y= '', tag = "C") +
  geom_point()

P4 <- ggplot(data_alphadiv, aes(x=ic, y=data_shannon)) +
  geom_boxplot(fill=c("blue","red")) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  geom_point()
P1
(P1 | P2) / (P3 | P4)

#### Por site, sintom?tico e controlo
P1 <- ggplot(data_alphadiv, aes(x=siteic, y=S.obs)) +
  geom_boxplot(fill=c("blue","red",'yellow','green')) +
  labs(title= 'Richness', x= ' ', y= '', tag = "A") +
  geom_point()

P2 <- ggplot(data_alphadiv, aes(x=siteic, y=S.chao1)) +
  geom_boxplot(fill=c("blue","red",'yellow','green')) +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  geom_point()

P3 <- ggplot(data_alphadiv, aes(x=siteic, y=data_evenness)) +
  geom_boxplot(fill=c("blue","red",'yellow','green')) +
  labs(title= 'Eveness', x= ' ', y= '', tag = "C") +
  geom_point()

P4 <- ggplot(data_alphadiv, aes(x=siteic, y=data_shannon)) +
  geom_boxplot(fill=c("blue","red",'yellow','green')) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  geom_point()
P1
(P1 | P2) / (P3 | P4)




#by tabanca
P1 <- ggplot(data_alphadiv, aes(x=Tabanca, y=S.obs)) +
  geom_boxplot(fill=c("blue","red", "black","pink", "yellow", 'green','white','purple','grey','orange')) +
  labs(title= 'Richness', x= ' ', y= '', tag = "A") +
  geom_point()

P2 <- ggplot(data_alphadiv, aes(x=Tabanca, y=S.chao1)) +
  geom_boxplot(fill=c("blue","red", "black","pink", "yellow",'green','white','purple','grey','orange')) +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  geom_point()

P3 <- ggplot(data_alphadiv, aes(x=Tabanca, y=data_evenness)) +
  geom_boxplot(fill=c("blue","red", "black","pink", "yellow",'green','white','purple','grey','orange')) +
  labs(title= 'Eveness', x= ' ', y= '', tag = "C") +
  geom_point()

P4 <- ggplot(data_alphadiv, aes(x=Tabanca, y=data_shannon)) +
  geom_boxplot(fill=c("blue","red", "black","pink", "yellow",'green','white','purple','grey','orange')) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  geom_point()
P1
(P1 | P2) / (P3 | P4)

#Statistical analyses

summary(aov(data_shannon ~ Site, data = data_alphadiv))
summary(aov(data_evenness ~ Site, data = data_alphadiv))

library(ggplot2)
# The mtcars dataset:
data <- as.matrix(data_otu1)

# Default Heatmap
heatmap(data, scale="column")

#####################################################################################

                        #Normalization per sample 

library(phyloseq)
library(dplyr)

data_phylo_filt <- filter_taxa(data_phylo, function(x) sum(x >= 4) >= (0.1 * length(x)), TRUE)
print(data_phylo_filt)



otu_table_filtered <- otu_table(data_phylo_filt)
total_counts_per_sample <- colSums(otu_table_filtered)
print(total_counts_per_sample)


otu_table_filtered1 <- otu_table(data_phylo_filt)
total_counts_per_sample1 <- rowSums(otu_table_filtered1)
print(total_counts_per_sample1)

total_counts_df <- as.data.frame(total_counts_per_sample1)
# Export as CSV
write.csv(total_counts_df, "total_counts_per_sample.csv", row.names = TRUE)

# Or, export as a generic text file with a custom separator (e.g., tab)
write.table(total_counts_df, "total_counts_per_sample.txt", sep = "\t", row.names = TRUE, quote = FALSE)
getwd()


# Access the OTU Table from your filtered phyloseq object
otu_table_filtered <- otu_table(data_phylo_filt)

# Calculate the number of OTUs per sample (non-zero counts)
num_otus_per_sample <- rowSums(otu_table_filtered > 0)

# Print or view the results
print(num_otus_per_sample)

total_counts_df <- as.data.frame(num_otus_per_sample)
# Export as CSV
write.csv(total_counts_df, "otus_per_sample.csv", row.names = TRUE)

# Or, export as a generic text file with a custom separator (e.g., tab)
write.table(total_counts_df, "otus_per_sample.txt", sep = "\t", row.names = TRUE, quote = FALSE)
getwd()



#In the further analyses, we will compare counts between samples. 
#To do so, we will have first to normalise the filtered occurrence table by sample 
#in order to obtain the same library size for every samples.

#Rarefy the data
set.seed(1782) # set seed for analysis reproducibility
OTU_filt_rar = rarefy_even_depth(otu_table(data_phylo), rngseed = TRUE, replace = FALSE) # rarefy the raw data using Phyloseq package
data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar)) # create a separated file

data_phylo_filt_rar <- phyloseq(OTU_filt_rar, SAM, TAX) # create a phyloseq object
data_phylo_filt_rar

rowSums(data_otu_filt_rar)

#####################################################################################
                                #Beta-diversity

#Distances calculation
#Beta-diversity represents the difference between two ecosystems/samples
#Distances metrics are between 0 and 1: 
#0 means identical communities in both samples and 1 means different communities in both samples.

# calculate Bray-Curtis distance using the vegan package
dist_bc <- as.matrix(vegdist(data_otu_filt_rar, method = "bray")) 

# a look at the first five rows / columns
dist_bc[1:5, 1:5]

#Visualisation using PCOA ordination plot

# calculate PCOA using Phyloseq package
pcoa_bc = ordinate(data_phylo_filt_rar, "PCoA", "bray") 

plot_ordination(data_phylo_filt_rar, pcoa_bc, color = "Site") + 
  geom_point(size = 3) + stat_ellipse(type='norm', linetype=2)

plot_ordination(data_phylo_filt_rar, pcoa_bc, color = "siteic") + 
  geom_point(size = 3) + stat_ellipse(type='norm', linetype=2)

plot_ordination(data_phylo_filt_rar, pcoa_bc, color = "Tabanca") + 
  geom_point(size = 3) + stat_ellipse(type='norm', linetype=2)

plot_ordination(data_phylo_filt_rar, pcoa_bc, color = "Tabanca") + 
  geom_point(size = 3)

plot_ordination(data_phylo_filt_rar, pcoa_bc, color = "Tabanca") + 
  geom_polygon(aes(fill=Tabanca)) + geom_point(size=5)

#The PCOA plot represents every samples as a dot, which is colored according to their sampling site.
#First, this two-dimensions PCOA plot show 52% of the total variance between the samples. 

