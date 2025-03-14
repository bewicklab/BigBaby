# Load necessary libraries
library(phyloseq)
library(microbiome)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library("tidyr")
library(tibble)
library(stringr)
#ibrary(ANCOMBC)

#Read biom file
data <- import_biom('ASV_Table.biom')

#Find OTU table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
tax_table_df<-tax_table(data)

metadata<-read.csv('Mole_MetadataAC4.csv',row.names=1, header= TRUE)

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#YOURFAMILYPHYLOSEQOBJECT<-tax_glom(YOURPHYLOSEQOBJECT,taxrank=‘Family’)

temp_phylo<-phyloseq(OTU_biom,tax_table_df)

#unique(tax_table(temp_phylo)[, "Rank5"]) #family
#temp_phylo<-tax_glom(temp_phylo,taxrank="Rank5")

Mphylo<-prune_samples(colnames(OTU_biom)[which(colnames(OTU_biom) %in% rownames(metadata))],temp_phylo)

M2phylo<-phyloseq(otu_table(Mphylo),tax_table(Mphylo),sample_data(metadata))

M2phylo = prune_taxa(assigned_taxa_seqs, M2phylo)

M2phylo = prune_taxa(taxa_names(M2phylo)[which(tax_table(M2phylo)[,1]=='k__Bacteria')],M2phylo)

# Find the minimum sequencing depth across all samples
min_depth <- min(sample_sums(M2phylo))

set.seed(1)

# Rarefy the entire dataset to the minimum depth
M2phylo <- rarefy_even_depth(M2phylo, sample.size = min_depth, rngseed = 123, verbose = TRUE)
sample_sums(M2phylo)

meta_data <- as.data.frame(sample_data(M2phylo))


### Prepare data for ANCOM; finding most abundant ASVs by class

#ggplot(meta_data, aes(x = StypeS, y = sample_sums(physeq_rarefied))) +
#  geom_boxplot() +
#  labs(x = "Class (StypeS)", y = "Sequencing Depth")

# Ensure sample IDs match
# Extract and orient OTU table properly
otu_table_df <- if (taxa_are_rows(M2phylo)) {
  as.data.frame(t(otu_table(M2phylo)))
} else {
  as.data.frame(otu_table(M2phylo))
}
seq_names <- colnames(otu_table_df)

row_sums <- rowSums(otu_table_df[names(otu_table_df) %in% seq_names])
sample_sums_values <- sample_sums(M2phylo)
all.equal(row_sums, sample_sums_values)  # Should return TRUE


# Add 'SampleID' as a column in the OTU table (using rownames)
otu_table_df$SampleID <- rownames(otu_table_df)

# Add `StypeS` from metadata to the OTU table
otu_table_df <- merge(otu_table_df, metadata["Gtype"], by.x = "SampleID", by.y = "row.names")

# Reshape OTU table into a long format for easier aggregation
asv_abundance_by_class <- otu_table_df %>%
  pivot_longer(cols = -c(SampleID, Gtype), names_to = "ASV", values_to = "Abundance") %>%
  group_by(Gtype, ASV) %>%
  summarize(TotalAbundance = sum(Abundance), .groups = "drop")


# Identify the top 100 ASVs for each `Gtype`
top_50_asvs_per_class <- asv_abundance_by_class %>%
  group_by(Gtype) %>%
  slice_max(order_by = TotalAbundance, n = 50)

# Convert to a list if desired
top_asvs_list <- split(top_50_asvs_per_class, top_50_asvs_per_class$Gtype)

### 50 most abundant ASVs overall...

total_abundance <- otu_table_df
row.names(total_abundance) <- total_abundance[,1]
total_abundance <- total_abundance[-1]
total_abundance <- data.frame(colSums(total_abundance[(names(total_abundance) != "Gtype")]))
total_abundance$ASV <- row.names(total_abundance)


#abundance_df <- rename(abundance_df, "TotalAbundance" = "colSums.total_abundance..names.total_abundance......Gtype....")
colnames(total_abundance)[1] <- "TotalAbundance"
#abundance_df$ASV <- row.names(abundance_df)

# Identify the 50 most abundant ASVs overall
top_50_asvs <- total_abundance %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 50)

# Print the top 50 ASVs
#print(top_50_asvs)

ASV_abun <- unique(c(top_asvs_list[["Larvae"]][["ASV"]], top_asvs_list[["Paedomorph"]][["ASV"]], top_asvs_list[["Metamorph"]][["ASV"]], top_50_asvs$ASV))

#---

# Ensure sample IDs match
row.names(otu_table_df) <- otu_table_df$SampleID
otu_table_df <- otu_table_df[row.names(meta_data), ]

#---

#shoving extra ASVs into its own ASV category
# Identify ASVs not in the keep_asvs vector
other_asvs <- setdiff(colnames(otu_table_df), ASV_abun)

# Create a new "Other" category by summing the reads for all ASVs not in keep_asvs
otu_table_df$other <- rowSums(otu_table_df[, other_asvs[!(other_asvs %in% c("SampleID", colnames(meta_data)))]])

# Subset the OTU table to keep only the specified ASVs and the new "Other" category
otu_table_reduced <- otu_table_df[, c(ASV_abun, "other")]



#---
#doesn't include 'other' taxon
otu_table_sub <- otu_table_df[names(otu_table_df) %in% ASV_abun] #missing 'other' ASVs if I want to run them


# Source the ANCOM script
source("D:/research/publications/facultative_mole_salamander_microbiomes/ANCOM/ancom2_w_tidyr.R")  # Replace with the path to the downloaded script

# Run ANCOM
ancom_res <-  ANCOM(t(otu_table_sub), meta_data, main_var="Gtype")

# View results
#print(ancom_res)

#save.image('ANCOM_ASV_moles.RData')

#--------------------------------------------------------

# Extract taxonomic information from phyloseq object
#tax_table_df <- as.data.frame(tax_table(M2phylo))  # Convert tax_table to a dataframe

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
tax_table_df<-tax_table(M2phylo)
colnames(tax_table_df)<-c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
tax_table_df<-gsub('p__', '', tax_table_df)
tax_table_df<-gsub('c__', '', tax_table_df)
tax_table_df<-gsub('o__', '', tax_table_df)
tax_table_df<-gsub('f__', '', tax_table_df)
tax_table_df<-gsub('g__', '', tax_table_df)
tax_table_df<-gsub('s__', '', tax_table_df)


for (k in 1:length(tax_table_df[,6])){
  if (is.na(tax_table_df[k,6])){
    if (!(is.na(tax_table_df[k,5]))){
      tax_table_df[k,6]<-paste0('unclassified_',tax_table_df[k,5])
    }else{
      if (!(is.na(tax_table_df[k,4]))){
        tax_table_df[k,6]<-paste0('unclassified_',tax_table_df[k,4])
      }else{
        if (!(is.na(tax_table_df[k,3]))){
          tax_table_df[k,6]<-paste0('unclassified_',tax_table_df[k,3])
        }
        else{
          tax_table_df[k,6]<-paste0('unclassified_',tax_table_df[k,2])
        }
      }
    }
  }
}
for (k in 1:length(tax_table_df[,6])){
  if (tax_table_df[k,6]=='NA'){
    if (tax_table_df[k,5]!='NA'){
      tax_table_df[k,6]<-paste0('unclassified_',tax_table_df[k,5])
    }else{
      if (tax_table_df[k,4]!='NA'){
        tax_table_df[k,6]<-paste0('unclassified_',tax_table_df[k,4])
      }else{
        if (tax_table_df[k,3]!='NA'){
          tax_table_df[k,6]<-paste0('unclassified_',tax_table_df[k,3])
        }
        else{
          tax_table_df[k,6]<-paste0('unclassified_',tax_table_df[k,2])
        }
      }
    }
  }
}

tax_table_df <- data.frame(tax_table_df)

tax_table_df <- tax_table_df %>%
  rownames_to_column(var = "taxa_id")  # Row names are the taxa_id values

#---

#remove?

# Replace "_NA" with "_Unclassified" in rank6
tax_table_df$Rank6 <- gsub("_NA$", "_Unclassified", tax_table_df$Rank6)

# Replace "_NA" with "_unclassified" in rank7
tax_table_df$Rank7 <- gsub("_NA$", "_unclassified", tax_table_df$Rank7)

tax_table_df <- tax_table_df %>%
  mutate(
    gen_sp = paste(
      str_remove(Rank6, "^g__"),  # Remove "g_" from Rank6
      str_remove(Rank7, "^s__"),  # Remove "s_" from Rank7
      sep = "_"
    )
  )

#---

tax_table_df$ASV <- paste(tax_table_df$gen_sp, tax_table_df$taxa_id, sep="_")


#---

# Step 1: Extract significant taxa based on W statistic
threshold <- 0.9 * max(ancom_res$out$W)  # Adjust threshold as needed

significant_taxa <- ancom_res$out %>%
  filter(W > threshold) %>% 
  select(taxa_id)

#extracting W-scores from ANCOM
w_scores <- ancom_res$out %>%
  select(taxa_id, W)  # Keep only taxa IDs and W-scores

# Merge taxonomic names with abundance data
w_scores <- w_scores %>%
  left_join(tax_table_df, by = "taxa_id")

#just keeping significant sequences
w_scores <- w_scores[w_scores$taxa_id %in% significant_taxa$taxa_id,] 
  ### THIS LINE MAYBE SHOULD BE REMOVED, I DONT THINK SO THOUGH!!!

# Order ASVs by descending W-score
ordered_taxa <- w_scores %>% arrange(desc(W)) %>% pull(taxa_id)

# Prepare data for W-score dot plot
wscore_data <- w_scores %>% 
  mutate(taxa_id = factor(taxa_id, levels = ordered_taxa))
#this w_scores object had 3 columns: "taxa_id", "W", "gen_sp"
### MAKE SURE THIS GETS SET UP CORRECTLY!!!

# Create a mapping of taxa_id to Rank7
wtaxa_to_asv <- wscore_data %>%
  distinct(taxa_id, ASV) %>%
  filter(taxa_id %in% ordered_taxa) %>%
  arrange(match(taxa_id, ordered_taxa)) # Ensure the order matches ordered_asvs

# Ensure that the taxa_id factor uses the same order as ordered_asvs
wscore_data <- wscore_data %>%
  mutate(
    taxa_id = factor(taxa_id, levels = ordered_taxa),
    asv_label = factor(taxa_id, levels = ordered_taxa, labels = wtaxa_to_asv$ASV))

#---

# Find the class with the highest abundance for each taxa_id
highest_abundance_class <- barplot_data %>%
  group_by(taxa_id) %>%
  slice_max(order_by = RelativeAbundance, n = 1, with_ties = FALSE) %>%  # Get the max RelativeAbundance per taxa_id
  select(taxa_id, Gtype) %>%
  rename(HighestAbundanceClass = Gtype)  # Rename for clarity

# Merge highest abundance class with W-scores
wscore_data <- wscore_data %>%
  left_join(highest_abundance_class, by = "taxa_id")

wscore_data <- wscore_data %>%
  mutate(HighestAbundanceClass = recode(HighestAbundanceClass, "Metamorph" = "Metamorphs"))

wscore_data <- wscore_data %>%
  mutate(HighestAbundanceClass = recode(HighestAbundanceClass, "Paedomorph" = "Paedomorphs"))

# Ensure taxa_id order matches ordered_asvs
wscore_data <- wscore_data %>%
  mutate(taxa_id = factor(taxa_id, levels = ordered_taxa))

#---

#ideally don't read in wonky pasted gen_sp/taxa_id labels, but the code works!
# Plot W-scores
ggplot(wscore_data, aes(x = W, y = asv_label, colour = HighestAbundanceClass)) +
  geom_point(size=5) +
  #scale_x_continuous(name = "W-Score") +
  scale_x_continuous(
    name = "W-Score",
    breaks = seq(floor(min(wscore_data$W)), ceiling(max(wscore_data$W)), by = 1)  # Generate breaks at intervals of 1
  ) +
  scale_y_discrete(
    name = "ASV",
    labels = function(labels) gsub("_[^_]*$", "", labels) # Remove text starting from the underscore
  ) +
  scale_colour_manual(
    name = "Life Stage",  # Legend title
    values = c("Larvae" = "seagreen1", "Paedomorphs" = "lightsteelblue1", "Metamorphs" = "coral4")  # Custom colors
  ) +
  labs(y = NULL) + #remove y-axis title
theme_bw()

#----------------------

#EXTRACTING RELATIVE ABUNDANCE DATA FOR STACKED BARPLOT

# Prepare data for stacked bar plot
relative_abundance <- transform_sample_counts(M2phylo, function(x) x / sum(x))
#not scaled to a percent...

barplot_data <- psmelt(relative_abundance) %>% 
  filter(OTU %in% ordered_taxa) %>% 
  mutate(taxa_id = factor(OTU, levels = ordered_taxa))

# Merge taxonomic names with abundance data
barplot_data <- barplot_data %>%
  left_join(tax_table_df, by = "taxa_id")

#### HOPEFULLY THIS WORKDS STILL
# Create a mapping of taxa_id to Rank7
taxa_to_asv <- barplot_data %>%
  distinct(taxa_id, ASV) %>%
  filter(taxa_id %in% ordered_taxa) %>%
  arrange(match(taxa_id, ordered_taxa)) # Ensure the order matches ordered_asvs

# Ensure that the taxa_id factor uses the same order as ordered_asvs
barplot_data <- barplot_data %>%
  mutate(
    taxa_id = factor(taxa_id, levels = ordered_taxa),
    asv_label = factor(taxa_id, levels = ordered_taxa, labels = taxa_to_asv$ASV)
  )


# Aggregate abundances by Gtype and Microbe (taxa_id), using the median
abundance_plot <- barplot_data %>%
  group_by(taxa_id, Gtype, asv_label) %>%
  summarize(MedianAbundance = median(Abundance), .groups = "drop") %>%
  group_by(taxa_id) %>%
  mutate(RelativeAbundance = MedianAbundance / sum(MedianAbundance) * 100) %>%
  ungroup()

# Ensure taxa_id order matches ordered_asvs
abundance_plot <- abundance_plot %>%
  mutate(taxa_id = factor(taxa_id, levels = ordered_taxa))

abundance_plot <- abundance_plot %>%
  mutate(Gtype = recode(Gtype, "Metamorph" = "Metamorphs"))

abundance_plot <- abundance_plot %>%
  mutate(Gtype = recode(Gtype, "Paedomorph" = "Paedomorphs"))

abundance_plot <- abundance_plot %>%
  mutate(Gtype = factor(Gtype, levels = c("Metamorphs", "Paedomorphs", "Larvae")))

#-----

#setting up occurrence data

# Ensure all combinations of taxa_id, Gtype, and asv_label exist
occurrence_plot <- barplot_data %>%
  # Ensure complete combinations of taxa_id, Gtype, and asv_label
  complete(
    asv_label = unique(barplot_data$asv_label),
    Gtype = unique(barplot_data$Gtype),
    taxa_id = unique(barplot_data$taxa_id),
    fill = list(Abundance = 0)
  ) %>%
  # Count occurrences where Abundance > 0
  group_by(Gtype, asv_label) %>%
  summarize(
    Occurrences = sum(Abundance > 0),
    taxa_id = first(taxa_id),  # Ensure taxa_id is preserved properly
    .groups = "drop"
  ) %>%
  # Calculate total samples for each Gtype
  left_join(
    barplot_data %>%
      distinct(Gtype, Sample) %>%
      group_by(Gtype) %>%
      summarize(TotalSamples = n(), .groups = "drop"),
    by = "Gtype"
  ) %>%
  # Calculate relative occurrences
  mutate(RelativeOccurrences = (Occurrences / TotalSamples) * 100) %>%
  ungroup()

#add w-scores to occurrence_plot
occurrence_plot <- occurrence_plot %>% 
  left_join(wscore_data[c("W", "asv_label")], by = "asv_label")

# Ensure taxa_id order matches ordered_asvs
occurrence_plot <- occurrence_plot %>%
  mutate(taxa_id = factor(taxa_id, levels = ordered_taxa))

occurrence_plot <- occurrence_plot %>%
  mutate(Gtype = recode(Gtype, "Metamorph" = "Metamorphs"))

occurrence_plot <- occurrence_plot %>%
  mutate(Gtype = recode(Gtype, "Paedomorph" = "Paedomorphs"))

occurrence_plot <- occurrence_plot %>%
  mutate(Gtype = factor(Gtype, levels = c("Metamorphs", "Paedomorphs", "Larvae")))


#---

#test plots

#this should be run with the tax labels, but not yet incorporated...
# Create the stacked bar plot with Rank7 as y-axis labels
ggplot(abundance_plot, aes(x = RelativeAbundance
                         , y = asv_label, fill = Gtype)) +
  geom_bar(stat = 'identity', position="stack") +
  #scale_y_discrete(name = "ASV") +
  scale_y_discrete(position = "right",
                   #name = "ASV",
                   labels = function(labels) gsub("_[^_]*$", "", labels),
                   # Remove text starting from the underscore
  ) +
  scale_x_continuous(name = "Relative Abundance (%)") +
  scale_fill_manual(
    name = "Life Stage",  # Legend title
    values = c("Larvae" = "seagreen1", "Paedomorphs" = "lightsteelblue1", "Metamorphs" = "coral4")  # Custom colors
  ) +
  labs(y = NULL) + #remove y-axis title
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Occurrence
ggplot(occurrence_plot, aes(x = RelativeOccurrences
                           , y = asv_label, fill = Gtype)) +
  geom_bar(stat = 'identity', position="stack") +
  #scale_y_discrete(name = "ASV") +
  #scale_y_discrete(position = "right",
  #                 #name = "ASV",
  #                 labels = function(labels) gsub("_[^_]*$", "", labels),
  #                 # Remove text starting from the underscore
  #) +
  scale_x_continuous(name = "Relative Occurrence (%)",
                     labels = function(x) paste0(round(x / 3), "%")  # Dynamically rescale labels
                     ) +
  scale_fill_manual(
    name = "Life Stage",  # Legend title
    values = c("Larvae" = "seagreen1", "Paedomorphs" = "lightsteelblue1", "Metamorphs" = "coral4")  # Custom colors
  ) +
  labs(y = NULL) + #remove y-axis title
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#---------------------------------------------------------------------------------------------------
#plots in paper

library("gridExtra")

A <- ggplot(occurrence_plot, aes(x = RelativeOccurrences, y = asv_label, fill = Gtype)) +
  geom_bar(stat = 'identity', position="stack") +
  #scale_y_discrete(name = "ASV") +
  #scale_y_discrete(position = "right",
  #                 #name = "ASV",
  #                 labels = function(labels) gsub("_[^_]*$", "", labels),
  #                 # Remove text starting from the underscore
  #) +
  scale_y_discrete(
    name = "W-Score",  # Y-axis title
    labels = setNames(occurrence_plot$W, occurrence_plot$asv_label)  # Replace y-axis labels with column W
  ) +
 scale_x_continuous(name = "Relative Occurrence (%)",
                    labels = function(x) paste0(round(x / 3), "%")) +
  scale_fill_manual(
    name = "Life Stage",  # Legend title
    values = c("Larvae" = "seagreen1", "Paedomorphs" = "lightsteelblue1", "Metamorphs" = "coral4")  # Custom colors
  ) +
 labs(y = NULL) + #remove y-axis title
  #scale_y_discrete(labels = NULL, breaks=NULL) +
  theme_bw(30) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold")) 
  
  
        
B <- ggplot(abundance_plot, aes(x = RelativeAbundance, y = asv_label, fill = Gtype)) +
  geom_bar(stat = 'identity', position="stack") +
  #scale_y_discrete(name = "ASV") +
  scale_y_discrete(position = "right",
                   #name = "ASV",
                   labels = function(labels) gsub("_[^_]*$", "", labels),
                   # Remove text starting from the underscore
  ) +
  scale_x_continuous(name = "Relative Abundance (%)") +
  scale_fill_manual(
    name = "Life Stage",  # Legend title
    values = c("Larvae" = "seagreen1", "Paedomorphs" = "lightsteelblue1", "Metamorphs" = "coral4")  # Custom colors
  ) +
  labs(y = NULL) + #remove y-axis title
  theme_bw(30) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "bold"))



library("gridExtra")
library("grid")
library("cowplot")

plots <- list(A, B)

plot_layout <- rbind(c(1,2))

grid_layout <- plot_grid(A, B, nrow = 1, rel_widths = c(1, 2.4))

ggsave("test3.png", grid_layout, height = 30, width = 20)



