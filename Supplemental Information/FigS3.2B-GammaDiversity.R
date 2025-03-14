library('vegan')
library('phyloseq')
library('ape')
library('microshades')
library('speedyseq')
library('ggplot2')
library('cowplot')
library('forcats')
library('tidyverse')
library('iNEXT.3D')

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

Mphylo<-prune_samples(sample_data(rMphylo)$StypeS=='Metamorph',rMphylo)
wPphylo<-prune_samples(sample_data(rMphylo)$StypeS=='w_Paedomorph',rMphylo)
sPphylo<-prune_samples(sample_data(rMphylo)$StypeS=='s_Paedomorph',rMphylo)
Lphylo<-prune_samples(sample_data(rMphylo)$StypeS=='Larvae',rMphylo)

Stages<-list(sign(otu_table(Mphylo)),sign(otu_table(wPphylo)),sign(otu_table(sPphylo)),sign(otu_table(Lphylo)))
names(Stages)<-c('Metamorph','w_Paedomorph','s_Paedomorph','Larvae')

gtemp<-iNEXT3D(Stages,diversity='TD',q=0,datatype='incidence_raw')


ggiNEXT3D(gtemp,type = 1, facet.var = "Order.q")+scale_color_manual(values=c('seagreen1','coral4','lightsteelblue1','royalblue1'))+scale_fill_manual(values=c('seagreen1','coral4','lightsteelblue1','royalblue1'))
