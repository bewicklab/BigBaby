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
library('StereoMorph')
library('tidyverse')
library('ggpubr')
library('MASS')
library('rlist')

#Function definition
Find_Projections<-function(Lc,Mc,Pc,Oc,Pp){
  orthoP<-orthogonalProjectionToLine(Pp, l1 = Lc, l2 = Mc)
  orthoP2<-orthogonalProjectionToLine(Pp, l1 = Pc, l2 = Oc)
  
  #Find the distances between the four points (larvae, metamorphs, paedomorph and paedomorph projection onto metamorphosis line)
  disttemp<-as.matrix(vegdist(as.matrix(rbind(Lc,Mc,Pp,orthoP)),method='euclidean',upper=TRUE,diag=TRUE))
  disttemp2<-as.matrix(vegdist(as.matrix(rbind(Oc,Pc,Pp,orthoP2)),method='euclidean',upper=TRUE,diag=TRUE))
  

  #The distance between the paedomorph projection and the larval centroid
  fromL_P<-disttemp[1,4]
  fromV_P<-disttemp2[1,4]
  #Scaled based on the distance between the larvae and metamorphs
  scaled_fromL_P<-disttemp[1,4]/disttemp[1,2]
  scaled_fromV_P<-disttemp2[1,4]/disttemp2[1,2]
  
  #The distance between the paedomorph projection and the metamorph centroid
  fromM_P<-disttemp[2,4]
  fromH_P<-disttemp2[2,4]
  #Scaled based on the distance between the larval and metamorph centroids
  scaled_fromM_P<-disttemp[2,4]/disttemp[1,2]
  scaled_fromH_P<-disttemp2[2,4]/disttemp2[1,2]
  
  #If we place the larvae at -0.5 and the metamorphs at 0.5, then...
  #At 0 the hybrid projection is exactly intermediate...
  #At negative values, the paedomorph projection is closer to the larvae
  #At positive values, the paedomorph projection is closer to the metamorph
  
  #If the hybrid projection lies between both parents...
  if (disttemp[1,4]/disttemp[1,2]<1 && disttemp[2,4]/disttemp[1,2]<1){
    metamorphosis_position_P<-0.5-disttemp[2,4]/disttemp[1,2]
  }else if (disttemp[1,4]>disttemp[2,4]){
    metamorphosis_position_P<-0.5+disttemp[2,4]/disttemp[1,2]
  }else{
    metamorphosis_position_P<--0.5-disttemp[1,4]/disttemp[1,2]
  }
  
  
  #If the hybrid projection lies between both parents...
  if (disttemp2[1,4]/disttemp2[1,2]<1 && disttemp2[2,4]/disttemp2[1,2]<1){
    metamorphosis_position_P2<-(1-disttemp2[2,4]/disttemp2[1,2])*1#disttemp2[1,2]/disttemp[1,2]
  }else if (disttemp2[1,4]>disttemp2[2,4]){
    metamorphosis_position_P2<-(1+disttemp2[2,4]/disttemp2[1,2])*1#disttemp2[1,2]/disttemp[1,2]
  }else{
    metamorphosis_position_P2<-(0-disttemp2[1,4]/disttemp2[1,2])*1#disttemp2[1,2]/disttemp[1,2]
  }
  
  
  output<-c(metamorphosis_position_P,metamorphosis_position_P2)
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


metadata<-read.csv('Mole_MetadataAC4.csv',row.names=1, header= TRUE)

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

temp_phylo<-phyloseq(OTU_biom,TAX,tree_file)
Mphylo<-prune_samples(colnames(OTU_biom)[which(colnames(OTU_biom) %in% rownames(metadata))],temp_phylo)

M2phylo<-phyloseq(otu_table(Mphylo),tax_table(Mphylo),sample_data(metadata),phy_tree(Mphylo))

M2phylo <- prune_taxa(assigned_taxa_seqs, M2phylo)

M2phylo <- prune_taxa(taxa_names(M2phylo)[which(tax_table(M2phylo)[,1]=='k__Bacteria')],M2phylo)

set.seed(1)

rMphylo <- rarefy_even_depth(M2phylo, rngseed=1, replace=F, verbose = TRUE)

#otu_table(rMphylo)<-sign(otu_table(rMphylo))

Mphylo<-prune_samples(sample_data(rMphylo)$StypeS=='Metamorph',rMphylo)
Lphylo<-prune_samples(sample_data(rMphylo)$StypeS=='Larvae',rMphylo)
#Using all larvae...
Lphylosmall<-Lphylo#prune_samples(sample_data(Lphylo)$svl<20,Lphylo)
Lphylobig<-Lphylo#prune_samples(sample_data(Lphylo)$svl>20,Lphylo)
Pphylo<-prune_samples(sample_data(rMphylo)$StypeS %in% c('w_Paedomorph','s_Paedomorph'),rMphylo)
Pphylobig<-prune_samples(sample_data(Pphylo)$svl>51,Pphylo)
Pphylosmall<-prune_samples(sample_data(Pphylo)$svl<52,Pphylo)
Pwphylo<-prune_samples(sample_data(rMphylo)$StypeS %in% c('w_Paedomorph'),rMphylo)
Psphylo<-prune_samples(sample_data(rMphylo)$StypeS %in% c('s_Paedomorph'),rMphylo)




Lcentroid<-as.vector(rowMeans(otu_table(Lphylosmall)))
Mcentroid<-as.vector(rowMeans(otu_table(Mphylo)))
Pcentroid<-as.vector(rowMeans(otu_table(Pphylobig)))
Ocentroid<-orthogonalProjectionToLine(Pcentroid, l1 = Lcentroid, l2 = Mcentroid)
Aphylo<-merge_phyloseq(Lphylo,Mphylo,Pphylosmall)
wAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='w_Paedomorph',Aphylo)
sAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='s_Paedomorph',Aphylo)
lAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='Larvae',Aphylo)
mAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='Metamorph',Aphylo)


#Dominant Taxa on Metamorphosis axis
temp<-Mcentroid-Lcentroid
summer<-0
for (k in 1:length(temp)){
  summer<-summer+temp[k]*temp[k]
}
summer<-sqrt(summer)
normalized_ancestral<-temp/summer
totalHeight<-sum(normalized_ancestral[which(normalized_ancestral>0)]*normalized_ancestral[which(normalized_ancestral>0)])
totalLine<-sum(normalized_ancestral[which(normalized_ancestral<0)]*normalized_ancestral[which(normalized_ancestral<0)])

metline_delta<-temp
metline_delta[metline_delta<=0]<-0
#Dominant Taxa on Paedomorph axis
metheight_delta<-temp
metheight_delta[metheight_delta>=0]<-0
metheight_delta_abs<-abs(metheight_delta)

#metline_delta_abs<-abs(metline_delta)
#temp<-sort(metline_delta_abs, index.return=TRUE,decreasing=TRUE)
#sorted_indices<-temp$ix
#sorted_metline_delta<-metline_delta[sorted_indices]
sorted_metline_delta<-metline_delta
sorted_metheight_delta<-metheight_delta

summer<-0
for (k in 1:length(sorted_metline_delta)){
  summer<-summer+sorted_metline_delta[k]*sorted_metline_delta[k]
}
summer<-sqrt(summer)

normalized_sorted_metline_delta<-sorted_metline_delta/summer

#temp<-sort(metheight_delta_abs, index.return=TRUE,decreasing=TRUE)
#sorted_indices<-temp$ix
#sorted_metheight_delta<-metheight_delta[sorted_indices]

summer<-0
for (k in 1:length(sorted_metheight_delta)){
  summer<-summer+sorted_metheight_delta[k]*sorted_metheight_delta[k]
}
summer<-sqrt(summer)

normalized_sorted_metheight_delta<-sorted_metheight_delta/summer

Ocentroid<-orthogonalProjectionToLine(sorted_metheight_delta, l1 = Lcentroid*0, l2 = sorted_metline_delta)

Ametline<-c()
Ametheight<-c()
for (k in 1:length(sample_names(Aphylo))){
  temp<-Find_Projections(Lcentroid*0,sorted_metline_delta,sorted_metheight_delta,Ocentroid,as.vector(otu_table(Aphylo)[,k])-Lcentroid)
  Ametline<-c(Ametline,temp[1])
  Ametheight<-c(Ametheight,temp[2])
}



df<-data.frame(sample_data(Aphylo)$StypeS,Ametline+0.5,Ametheight)
colnames(df)<-c('stype','metline','metheight')

par(mar=c(5,5,5,2))
#plot(df$metline,df$metheight,xlab='Ancestral (Metamorph) Axis',ylab='Loss of Larval Taxa',cex.lab=1.5)
plot(c(0,1),c(0,1),xlab=paste('Gain of Metamorph Taxa',paste0('(',round(totalLine,3)*100,'%)')),ylab=paste('Loss of Larval Taxa',paste0('(',round(totalHeight,3)*100,'%)')),cex.lab=1.5,ylim=c(-0.25,1.1),xlim=c(-0.25,1.1),cex=2)
#points(df$metline[which(df$stype=='Metamorph')],df$metheight[which(df$stype=='Metamorph')],pch=16,col='coral4')
points(df$metline[which(df$stype=='w_Paedomorph')],df$metheight[which(df$stype=='w_Paedomorph')],pch=16,col='royalblue1')
points(df$metline[which(df$stype=='s_Paedomorph')],df$metheight[which(df$stype=='s_Paedomorph')],pch=16,col='lightsteelblue1')
lines(c(-0.5,3),c(0,0),lty='dashed',col='lightgrey')
lines(c(0,0),c(-1,2),lty='dashed',col='lightgrey')
lines(c(-0.5,3),c(-0.5,3))
lines(c(0.5,-1),c(0.5,2))
points(0,0,pch=16,col='seagreen1',cex=2)
points(1,1,pch=16,col='coral4',cex=2)
points(c(0,1),c(0,1),cex=2)
pp1<-c(mean(df$metline[which(df$stype=='w_Paedomorph')]),mean(df$metheight[which(df$stype=='w_Paedomorph')]))
ln1<-orthogonalProjectionToLine(pp1, l1 = c(0,0), l2 = c(1,1))
lines(c(pp1[1],ln1[1]),c(pp1[2],ln1[2]),lty='dashed',col='lightgrey')
pp2<-c(mean(df$metline[which(df$stype=='s_Paedomorph')]),mean(df$metheight[which(df$stype=='s_Paedomorph')]))
ln2<-orthogonalProjectionToLine(pp2, l1 = c(0,0), l2 = c(1,1))
lines(c(pp2[1],ln2[1]),c(pp2[2],ln2[2]),lty='dashed',col='lightgrey')

ln3<-orthogonalProjectionToLine(pp1, l1 = c(0.5,0.5), l2 = c(0,1))
lines(c(pp1[1],ln3[1]),c(pp1[2],ln3[2]),lty='dashed',col='lightgrey')
ln4<-orthogonalProjectionToLine(pp2, l1 = c(0.5,0.5), l2 = c(0,1))
lines(c(pp2[1],ln4[1]),c(pp2[2],ln4[2]),lty='dashed',col='lightgrey')

points(mean(df$metline[which(df$stype=='w_Paedomorph')]),mean(df$metheight[which(df$stype=='w_Paedomorph')]),pch=16,col='royalblue1',cex=2)
points(mean(df$metline[which(df$stype=='s_Paedomorph')]),mean(df$metheight[which(df$stype=='s_Paedomorph')]),pch=16,col='lightsteelblue1',cex=2)
points(mean(df$metline[which(df$stype=='w_Paedomorph')]),mean(df$metheight[which(df$stype=='w_Paedomorph')]),cex=2)
points(mean(df$metline[which(df$stype=='s_Paedomorph')]),mean(df$metheight[which(df$stype=='s_Paedomorph')]),cex=2)




