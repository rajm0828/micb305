set.seed(123)
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(picante)

#Load Data 
metafp <- "data/hiv_covid_metadata_fix.tsv"
meta <- read_delim(metafp, delim="\t")
otufp <- "data/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)
taxfp <- "data/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")
phylotreefp <- "data/tree.nwk"
phylotree <- read.tree(phylotreefp)

#Format OTU table
# OTU tables should be a matrix with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back
#Save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
#Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
#Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

# Format sample metadata
#Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
#Make sample ids the rownames
rownames(samp_df)<- meta$'sample-id'
#Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
view(SAMP)

#Format taxonomy
#Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
#Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
#Make sample ids the rownames
rownames(tax_mat) <- tax$`Feature ID`
#Make taxa table
TAX <- tax_table(tax_mat)

class(OTU)
class(SAMP)
class(TAX)
view(TAX)

#Create phyloseq object: Merge all into a phyloseq object
hiv_covid <- phyloseq(OTU, SAMP, TAX, phylotree)

#Filter the phyloseq to remove chloroplast and mitochondria samples
hiv_bacteria_only <- subset_taxa(hiv_covid, Domain == "d__Bacteria" & Class != "c__Chloroplast"
                                 & Family != "f__Mitochondria")

#Filter HIV Cohort
hiv_bacteria_a <- subset_samples(hiv_bacteria_only, Cohort_Short == "A" | is.na(Cohort_Short))

#Filter HIV Visit
hiv_bacteria_a_2 <- subset_samples(hiv_bacteria_a, Visit == 2 | is.na(Visit))

#Filter COVID Stool
hiv_bacteria_a_2_stool <- subset_samples(hiv_bacteria_a_2, env_medium == "Stool" | is.na(env_medium))

#Filter COVID + ICU
hiv_covid_filter_final <- subset_samples(hiv_bacteria_a_2_stool, Host_disease %in% c("Covid", "ICU") | is.na(Host_disease))

#View Final Filtered Data
sample_data(hiv_covid_filter_final)
ntaxa(hiv_covid_filter_final)

#Generate alpha rarefaction curve
hiv_covid_curve <- rarecurve(t(as.data.frame(otu_table(hiv_covid_filter_final))), cex=0.1)

#Generate rarefied phyloseq object
hiv_covid_rare <- rarefy_even_depth(hiv_covid_filter_final, rngseed = 1, sample.size = 45366)

#Generate Alpha Diversity Metrics
gg_richness <- plot_richness(hiv_covid_rare, x = "origin") +
  xlab("Subject ID") +
  geom_boxplot()
gg_richness

#Generate Phylogenetic Diversity

#Calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(hiv_covid_rare)), phy_tree(hiv_covid_rare),
                 include.root=F) 

#Add PD to metadata table
sample_data(hiv_covid_rare)$PD <- phylo_dist$PD

#Plot any metadata category against the PD
plot.pd <- ggplot(sample_data(hiv_covid_rare), aes(origin, PD)) + 
  geom_boxplot() +
  xlab("Subject ID") +
  ylab("Phylogenetic Diversity")

#View Phylogenetic Diversity Plot
plot.pd

#Generate Beta Diversity Metrics
#Bray
bc_dm <- distance(hiv_covid_rare, method="bray")

pcoa_bc <- ordinate(hiv_covid_rare, method="PCoA", distance=bc_dm)

plot_ordination(hiv_covid_rare, pcoa_bc, color = "origin")

#Jaccard
jc_dm <- distance(hiv_covid_rare, method="jaccard")

pcoa_jc <- ordinate(hiv_covid_rare, method="PCoA", distance=jc_dm)

plot_ordination(hiv_covid_rare, pcoa_jc, color = "origin")

#Unweighted Unifrac
un_dm <- distance(hiv_covid_rare, method="uunifrac")

pcoa_un <- ordinate(hiv_covid_rare, method="PCoA", distance=un_dm)

plot_ordination(hiv_covid_rare, pcoa_un, color = "origin")

#Weighted Unifrac
wu_dm <- distance(hiv_covid_rare, method="wunifrac")

pcoa_wu <- ordinate(hiv_covid_rare, method="PCoA", distance=wu_dm)

plot_ordination(hiv_covid_rare, pcoa_wu, color = "origin")


