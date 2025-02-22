library('vegan')
library('phyloseq')
library('ape')
library('microshades')
library('speedyseq')
library('ggplot2')
library('cowplot')
library('forcats')
library('tidyverse')

#Modified function for creating the custom legend
custom_legend2column <- function (mdf, cdf, group_level = "Phylum", subgroup_level = "Genus", x = "Sample",
                                  y = "Abundance", legend_key_size = 0.4, legend_text_size = 10)
{
  if (is.null(mdf[[group_level]])) {
    stop("mdf 'group_level' does not exist")
  }
  
  if (is.null(mdf[[subgroup_level]])) {
    stop("mdf 'subgroup_level' does not exist")
  }
  
  if (is.null(cdf$hex)) {
    stop("cdf 'hex' does not exist")
  }
  
  col_name_group <- paste0("Top_", group_level)
  col_name_subgroup <- paste0("Top_", subgroup_level)
  
  group_level_names <- unique(cdf[[col_name_group]])
  
  for (i in 1:length(group_level_names))
  {
    if( i == 1)
    {
      complete_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size = legend_text_size)
      tracer<-length(unique(mdf[which(mdf[,21]==group_level_names[i]),22]))+2
    }
    else
    {
      new_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size =legend_text_size)
      
      complete_height <- tracer
      new_height <-length(unique(mdf[which(mdf[,21]==group_level_names[i]),22]))+2
      tracer<-tracer+new_height
      print(c(complete_height,new_height))
      
      complete_legend <-plot_grid(complete_legend, new_legend, ncol=1, rel_heights = c(complete_height,new_height))
    }
  }
  plot(complete_legend)
  complete_legend
}

individual_legend2 <- function (mdf,
                                cdf,
                                group_name,
                                col_name_group = "Top_Phylum",
                                col_name_subgroup = "Top_Genus",
                                x = "Sample",
                                y = "Abundance",
                                legend_key_size = 0.4,
                                legend_text_size = 10)
{
  select_mdf <- mdf %>% filter(!!sym(col_name_group) == group_name)
  select_cdf <- cdf %>% filter(!!sym(col_name_group) == group_name)
  
  select_plot <- ggplot(select_mdf,
                        aes_string(x = x, y = y, fill = col_name_subgroup, text = col_name_subgroup)) +
    geom_col( position="fill") +
    scale_fill_manual(name = group_name,
                      values = select_cdf$hex,
                      breaks = select_cdf[[col_name_subgroup]]) +
    theme(legend.justification = "left") +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.key.size = unit(legend_key_size, "lines"), text=element_text(size=legend_text_size))
  
  legend <- get_legend(select_plot)
}


#Read biom file
data <- import_biom('ASV_Table.biom')

#Find OTU table in biom file
OTU_biom<-otu_table(data)


#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-multi2di(read.tree('Drtree.nwk'))

#Rename the tree tips to match the ASV table names
cut_names<-c()
for (k in 1:length(taxa_names(tree_file))){
  #Depending on the Qiime2 output, you may need to split the names with a space or with an underscore
  #cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],' ')[[1]][1])
  cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],'_')[[1]][1])
}
taxa_names(tree_file)<-cut_names

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)
colnames(TAX)<-c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
TAX<-gsub('p__', '', TAX)
TAX<-gsub('c__', '', TAX)
TAX<-gsub('o__', '', TAX)
TAX<-gsub('f__', '', TAX)
TAX<-gsub('g__', '', TAX)
TAX<-gsub('s__', '', TAX)


for (k in 1:length(TAX[,6])){
  if (is.na(TAX[k,6])){
    if (!(is.na(TAX[k,5]))){
      TAX[k,6]<-paste0('unclassified_',TAX[k,5])
    }else{
      if (!(is.na(TAX[k,4]))){
        TAX[k,6]<-paste0('unclassified_',TAX[k,4])
      }else{
        if (!(is.na(TAX[k,3]))){
          TAX[k,6]<-paste0('unclassified_',TAX[k,3])
        }
        else{
          TAX[k,6]<-paste0('unclassified_',TAX[k,2])
        }
      }
    }
  }
}
for (k in 1:length(TAX[,6])){
  if (TAX[k,6]=='NA'){
    if (TAX[k,5]!='NA'){
      TAX[k,6]<-paste0('unclassified_',TAX[k,5])
    }else{
      if (TAX[k,4]!='NA'){
        TAX[k,6]<-paste0('unclassified_',TAX[k,4])
      }else{
        if (TAX[k,3]!='NA'){
          TAX[k,6]<-paste0('unclassified_',TAX[k,3])
        }
        else{
          TAX[k,6]<-paste0('unclassified_',TAX[k,2])
        }
      }
    }
  }
}


metadata<-read.csv('Mole_MetadataAC4.csv',row.names=1, header= TRUE)

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]


temp_phylo<-phyloseq(OTU_biom,TAX,tree_file)
Mphylo<-prune_samples(colnames(OTU_biom)[which(colnames(OTU_biom) %in% rownames(metadata))],temp_phylo)

M2phylo<-phyloseq(otu_table(Mphylo),tax_table(Mphylo),sample_data(metadata),phy_tree(Mphylo))

M2phylo = prune_taxa(assigned_taxa_seqs, M2phylo)

M2phylo = prune_taxa(taxa_names(M2phylo)[which(tax_table(M2phylo)[,1]=='k__Bacteria')],M2phylo)

set.seed(1)

rMphylo = rarefy_even_depth(M2phylo, rngseed=1, replace=F, verbose = TRUE)
sample_data(rMphylo)$StypeS[which(sample_data(rMphylo)$StypeS=='Larvae')]<-'Larvae'
sample_data(rMphylo)$StypeS[which(sample_data(rMphylo)$StypeS=='s_Paedomorph')]<-'s_Paedomorphs'
sample_data(rMphylo)$StypeS[which(sample_data(rMphylo)$StypeS=='w_Paedomorph')]<-'w_Paedomorphs'
sample_data(rMphylo)$StypeS[which(sample_data(rMphylo)$StypeS=='Metamorph')]<-'Metamorphs'

#%%%%%%%%%%%%%%%%%%%%%%   MAKE FANCY BAR GRAPHS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum
summed_phyla<-rowSums(data.frame(otu_table(tax_glom(rMphylo,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list<-data.frame(tax_table(tax_glom(rMphylo,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_phyla_list<-phyla_list[order(-summed_phyla)]

#Prepare the OTU table in the correct format for the plotting program (this must be done with level 6 data including genera!)
mdf_prep <- prep_mdf(rMphylo)
#Tell the function to plot the five most abundant phyla (you could pick something different... )
color_objs_GP <- create_color_dfs(mdf_prep,selected_groups = c('Proteobacteria','Bacteroidetes','Actinobacteria','Firmicutes','Gracilibacteria'),cvd=TRUE)
#Extract the OTU table and color choices (again, putting things in the right format for the function)
mdf_GP <- color_objs_GP$mdf
cdf_GP <- color_objs_GP$cdf
new_cdf_GP <- color_reassign(cdf_GP,
                          group_assignment = c("Bacteroidetes", "Gracilibacteria","Firmicutes","Proteobacteria","Actinobacteria"),
                          color_assignment = c("micro_orange", "micro_brown","micro_purple","micro_cvd_green","micro_cvd_blue"))

#Expand the number of Actinobacteria shown
new_groups <- extend_group(mdf_GP, new_cdf_GP, "Phylum", "Genus", "Proteobacteria", existing_palette = "micro_cvd_green", new_palette = "micro_green", n_add = 5)

#Define your new legend with the expanded groups (FUNCTIONS TO GENERATE THIS CUSTOM LEGEND ARE AT THE BOTTOM OF THIS R FILE)
GP_legend_new <-custom_legend2column(new_groups$mdf, new_groups$cdf,legend_key_size=0.6,legend_text_size = 10)

#Define the plot that you will be making
plot2 <- plot_microshades(new_groups$mdf, new_groups$cdf)
#Define all the formatting aspects of the plot
plot_diff2 <- plot2 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 6)) +
  facet_grid(~fct_relevel(StypeS,'Metamorphs','w_Paedomorphs','s_Paedomorphs','Larvae'), scale="free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.text.x = element_text(size= 4)) +
  theme(plot.margin = margin(6,20,6,6))+
  theme(axis.title.x = element_text(size= 20))+
  theme(axis.title.y = element_text(size= 20))

#Make the plot
plot_grid(plot_diff2, GP_legend_new,  rel_widths = c(1, .25))


#%%%%%%%%%%%%%%%%%%%%%%   FIND PHYLUM PERCENTAGES ACROSS ALL SPECIES/POPULATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_phyla_order<-summed_phyla[order(-summed_phyla)]
#Turn the vector of sums into a dataframe
abundant_phyla<-data.frame(summed_phyla_order[1:5]*100/sum(summed_phyla))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_phyla)<-sorted_phyla_list[1:5]



#%%%%%%%%%%%%%%%%%%%%%%   FUNCTIONS FOR MAKING NICELY SCALED LEGENDS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


custom_legend2column <- function (mdf, cdf, group_level = "Phylum", subgroup_level = "Genus", x = "Sample",
                                  y = "Abundance", legend_key_size = 0.4, legend_text_size = 10)
{
  if (is.null(mdf[[group_level]])) {
    stop("mdf 'group_level' does not exist")
  }
  
  if (is.null(mdf[[subgroup_level]])) {
    stop("mdf 'subgroup_level' does not exist")
  }
  
  if (is.null(cdf$hex)) {
    stop("cdf 'hex' does not exist")
  }
  
  col_name_group <- paste0("Top_", group_level)
  col_name_subgroup <- paste0("Top_", subgroup_level)
  
  group_level_names <- unique(cdf[[col_name_group]])
  
  for (i in 1:length(group_level_names))
  {
    if( i == 1)
    {
      complete_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size = legend_text_size)
      tracer<-length(unique(mdf[which(mdf[,21]==group_level_names[i]),22]))+2
    }
    else
    {
      new_legend <-individual_legend2 (mdf, cdf, group_level_names[i], col_name_group, col_name_subgroup, legend_key_size = legend_key_size, legend_text_size =legend_text_size)
      
      complete_height <- tracer
      new_height <-length(unique(mdf[which(mdf[,21]==group_level_names[i]),22]))+2
      tracer<-tracer+new_height
      
      complete_legend <-plot_grid(complete_legend, new_legend, ncol=1, rel_heights = c(complete_height,new_height))
    }
  }
  complete_legend
}

individual_legend2 <- function (mdf,
                                cdf,
                                group_name,
                                col_name_group = "Top_Phylum",
                                col_name_subgroup = "Top_Genus",
                                x = "Sample",
                                y = "Abundance",
                                legend_key_size = 0.4,
                                legend_text_size = 10)
{
  select_mdf <- mdf %>% filter(!!sym(col_name_group) == group_name)
  select_cdf <- cdf %>% filter(!!sym(col_name_group) == group_name)
  
  select_plot <- ggplot(select_mdf,
                        aes_string(x = x, y = y, fill = col_name_subgroup, text = col_name_subgroup)) +
    geom_col( position="fill") +
    scale_fill_manual(name = group_name,
                      values = select_cdf$hex,
                      breaks = select_cdf[[col_name_subgroup]]) +
    theme(legend.justification = "left") +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.key.size = unit(legend_key_size, "lines"), text=element_text(size=legend_text_size))
  
  legend <- get_legend(select_plot)
}

