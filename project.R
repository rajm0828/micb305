#Set the seed to make program reproducible
set.seed(123)

#Load libraries
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(picante)
library(DESeq2)
library(ggsignif)
library(indicspecies)

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


#### Alpha and beta diversity ####

#Extract diversity metrics
alpha_df <- estimate_richness(hiv_covid_rare, measures = c("Observed", "Shannon"))

#Extract metadata from the rarefied phyloseq object and combine
meta_rare <- data.frame(sample_data(hiv_covid_rare))
meta_rare <- cbind(meta_rare, alpha_df)

#Generate Alpha Diversity Metrics - Observed
gg_observed <- plot_richness(hiv_covid_rare, x = "origin", measures = "Observed") +
  xlab("Subject ID") +
  geom_boxplot() +
  geom_signif(comparisons = list(c("covid_dataset", "hiv_dataset")), annotation = "*") +
  theme_classic(base_size = 12) 
gg_observed

#Find Observed Significance
hist(meta_rare$Observed)
shapiro.test(meta_rare$Observed)
#Observed is non-parametric
wilcox.test(Observed ~ origin, data = meta_rare, exact = FALSE)

#Generate Alpha Diversity Metrics - Shannon
gg_shannon <- plot_richness(hiv_covid_rare, x = "origin", measures = "Shannon") +
  xlab("Subject ID") +
  geom_boxplot() +
  geom_signif(comparisons = list(c("covid_dataset", "hiv_dataset")), annotation = "*") +
  theme_classic(base_size = 12) 
gg_shannon

#Find Shannon Significance
hist(meta_rare$Shannon)
shapiro.test(meta_rare$Shannon)
#Shannon is non-parametric
wilcox.test(Shannon ~ origin, data = meta_rare, exact = FALSE)

#Generate Beta Diversity Metrics
#Bray
bc_dm <- phyloseq::distance(hiv_covid_rare, method = "bray")
pcoa_bc <- ordinate(hiv_covid_rare, method="PCoA", distance=bc_dm)
plot_ordination(hiv_covid_rare, pcoa_bc, color = "origin") +
  stat_ellipse() +
  scale_color_manual(values = c("hiv_dataset" = "#D97979", "covid_dataset" = "#6495C8")) +

  theme_classic(base_size = 12)

#Find Bray Significance
dm_braycurtis <- vegdist(t(otu_table(hiv_covid_rare)), method="bray") # Bray-curtis
dat <- as.data.frame(sample_data(hiv_covid_rare))
adonis2(dm_braycurtis ~ origin, data = meta_rare)


#### Taxonomy bar plots ####

# Convert to relative abundance
hiv_covid_rare_RA <- transform_sample_counts(hiv_covid_rare, function(x) x/sum(x))

# To remove black bars, "glom"/group by phylum first
hiv_covid_phylum <- tax_glom(hiv_covid_rare_RA, taxrank = "Phylum", NArm=FALSE)

hiv_covid_taxa_plot <- plot_bar(hiv_covid_phylum, fill="Phylum") + 
  facet_wrap(.~origin, scales = "free_x")
hiv_covid_taxa_plot


#### Differential abundance ####

# converting phyloseq object into deseq
phyloseq_object_plus1 <- transform_sample_counts(hiv_covid_filter_final, function(x) x+1)
deseq_phylo <- phyloseq_to_deseq2(phyloseq_object_plus1, ~origin)
deseq_object <- DESeq(deseq_phylo)

results <- results(deseq_object, tidy=TRUE)
view(results)

# visualizing as a volcano plot
ggplot(results) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

vol_plot <- results |>
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot

# visualizing as a bar plot

# new table of results for bar plot
sigASVs <- results |> 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)

# ASV names
sigASVs_vec <- sigASVs|>
  pull(ASV)

# prune 
pruned_deseq <- prune_taxa(sigASVs_vec,hiv_covid_filter_final)
sigASVs <- tax_table(TAX) |>
  as.data.frame() |>
  rownames_to_column(var="ASV") |>
  right_join(sigASVs) |>
  arrange(log2FoldChange) |>
  mutate(Family = make.unique(Family)) |>
  mutate(Family = factor(Family, levels=unique(Family)))

# barplot
ggplot(sigASVs) +
  geom_bar(aes(x=Family, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Family, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))


#### Indicator taxa analysis ####

# glom to Genus
hiv_covid_genus <- tax_glom(hiv_covid_filter_final, "Genus", NArm = FALSE)

# convert to relative abundance
hiv_covid_genus_RA <- transform_sample_counts(hiv_covid_genus, fun=function(x) x/sum(x))

#ISA
isa_hiv_covid <- multipatt(t(otu_table(hiv_covid_genus_RA)), cluster = sample_data(hiv_covid_genus_RA)$`origin`)
summary(isa_hiv_covid)

## visualize as a table with ASV names and taxonomy ##
isa_taxtable <- tax_table(hiv_covid_filter_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_hiv_covid_table <- isa_hiv_covid$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(isa_taxtable) %>%
  filter(p.value<=0.05)
View(isa_hiv_covid_table)

write.csv(isa_hiv_covid_table, "isa_table_full.csv", row.names = FALSE)

## visualize as a bar plot ##
# prepping data
isa_hiv_covid_plot <- isa_hiv_covid$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(isa_taxtable, by = "ASV") %>%
  filter(p.value <= 0.05) %>%
  mutate(
    Group = case_when(
      s.hiv_dataset == 1 ~ "HIV",
      s.covid_dataset == 1 ~ "COVID"
    ),
    Taxon = if_else(
      !is.na(Genus) & Genus != "",
      Genus,
      ASV
    )
  )

# selecting the top 20 indicators per group
isa_hiv_covid_top20 <- isa_hiv_covid_plot %>%
  group_by(Group) %>%
  slice_max(order_by = stat, n = 20) %>%
  ungroup()

write.csv(isa_hiv_covid_top20, "isa_table_top20.csv", row.names = FALSE)

# making the bar plot
isa_hiv_covid_plot_final <- ggplot(isa_hiv_covid_top20,
       aes(x = stat, y = reorder(Taxon, stat), fill = Group)) +
  geom_col() +
  facet_wrap(~Group, scales = "free_y") +
  scale_fill_manual(
    values = c("HIV" = "#D97979", "COVID" = "#6495C8")
  ) +
  labs(
    x = "Indicator Value",
    y = "Genus",
    fill = "Infection Group"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  )
isa_hiv_covid_plot_final
