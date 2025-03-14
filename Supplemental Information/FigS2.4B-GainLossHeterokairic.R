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
  s1<-sign(sum(orthoP*Mc))
  s2<-sign(sum(orthoP2*Oc))
  
  
  #Find the distances between the four points (larvae, metamorphs, paedomorph and paedomorph projection onto metamorphosis line)
  disttemp<-as.matrix(vegdist(as.matrix(rbind(Lc,Mc,Pp,orthoP)),method='euclidean',upper=TRUE,diag=TRUE))
  disttemp2<-as.matrix(vegdist(as.matrix(rbind(Pc,Oc,Pp,orthoP2)),method='euclidean',upper=TRUE,diag=TRUE))
  

  
 
  metamorphosis_position_P<-s1*disttemp[1,4]/disttemp2[1,2]
  
  #If the hybrid projection lies between both parents...
  #if (disttemp2[1,4]/disttemp2[1,2]<1 && disttemp2[2,4]/disttemp2[1,2]<1){
  #  metamorphosis_position_P2<-(1-disttemp2[2,4]/disttemp2[1,2])*1#disttemp2[1,2]/disttemp[1,2]
  #}else if (disttemp2[1,4]>disttemp2[2,4]){
  #  metamorphosis_position_P2<-(1+disttemp2[2,4]/disttemp2[1,2])*1#disttemp2[1,2]/disttemp[1,2]
  #}else{
  #  metamorphosis_position_P2<-(0-disttemp2[1,4]/disttemp2[1,2])*1#disttemp2[1,2]/disttemp[1,2]
  #}
  metamorphosis_position_P2<-s2*disttemp2[1,4]/disttemp2[1,2]
  
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


metadata<-read.csv('1METADATAa.csv',row.names=1, header= TRUE)

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


Mphylo<-prune_samples(sample_data(rMphylo)$StypeS=='Metamorph',rMphylo)
Lphylo<-prune_samples(sample_data(rMphylo)$StypeS=='Larvae',rMphylo)
#Using all larvae...
Lphylosmall<-Lphylo#prune_samples(sample_data(Lphylo)$svl<20,Lphylo)
Lphylobig<-Lphylo#prune_samples(sample_data(Lphylo)$svl>20,Lphylo)
Pphylo<-prune_samples(sample_data(rMphylo)$StypeS %in% c('Paedomorphwinter'),rMphylo)
Pphylobig<-prune_samples(sample_data(Pphylo)$svl>51,Pphylo)
Pphylosmall<-prune_samples(sample_data(Pphylo)$svl<52,Pphylo)
Pwphylo<-prune_samples(sample_data(rMphylo)$StypeS %in% c('Paedomorphwinter'),rMphylo)



Lcentroid<-as.vector(rowMeans(otu_table(Lphylosmall)))
Mcentroid<-as.vector(rowMeans(otu_table(Mphylo)))
Pcentroid<-as.vector(rowMeans(otu_table(Pphylobig)))
Ocentroid<-orthogonalProjectionToLine(Pcentroid, l1 = Lcentroid, l2 = Mcentroid)
Aphylo<-merge_phyloseq(Lphylo,Mphylo,Pphylosmall)
wAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='Paedomorphwinter',Aphylo)
lAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='Larvae',Aphylo)
mAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='Metamorph',Aphylo)


#Dominant Taxa on Metamorphosis axis
temp<-Pcentroid-Ocentroid
temp2<-temp
summer<-0
for (k in 1:length(temp)){
  summer<-summer+temp2[k]*temp2[k]
}
summer<-sqrt(summer)
temp2<-temp2/summer

totalLine<-sum(temp2[which(temp2<0)]*temp2[which(temp2<0)])
totalHeight<-sum(temp2[which(temp2>0)]*temp2[which(temp2>0)])

metline_delta<-temp
metline_delta[metline_delta>=0]<-0
#Dominant Taxa on Paedomorph axis
metheight_delta<-temp
metheight_delta[metheight_delta<=0]<-0
metheight_delta_abs<-abs(metheight_delta)

#metline_delta_abs<-abs(metline_delta)
#temp<-sort(metline_delta_abs, index.return=TRUE,decreasing=TRUE)
#sorted_indices<-temp$ix
#sorted_metline_delta<-metline_delta[sorted_indices]
sorted_metline_delta<-metline_delta
sorted_metheight_delta<-metheight_delta

Ametline<-c()
Ametheight<-c()
for (k in 1:length(sample_names(Aphylo))){
  temp<-Find_Projections(Lcentroid*0,sorted_metline_delta,Lcentroid*0,sorted_metheight_delta,as.vector(otu_table(Aphylo)[,k])-Ocentroid)
  Ametline<-c(Ametline,temp[1])
  Ametheight<-c(Ametheight,temp[2])
}



df<-data.frame(sample_data(Aphylo)$StypeS,Ametline,Ametheight)
colnames(df)<-c('stype','metline','metheight')

Lx<-mean(df$metline[which(df$stype=='Larvae')])
Ly<-mean(df$metheight[which(df$stype=='Larvae')])
Mx<-mean(df$metline[which(df$stype=='Metamorph')])
My<-mean(df$metheight[which(df$stype=='Metamorph')])
Px<-mean(df$metline[which(df$stype=='Paedomorphwinter')])
Py<-mean(df$metheight[which(df$stype=='Paedomorphwinter')])

ff<-Find_Projections(Lcentroid*0,sorted_metline_delta,Lcentroid*0,sorted_metheight_delta,Pcentroid-Ocentroid)

par(mar=c(5,5,5,2))
plot(c(-10,-10),c(-10,-11),xlab=paste('Loss of Ancestral Taxa',paste0('(',round(totalLine,3)*100,'.0%)')),ylab=paste('Gain of Heterokairic Taxa',paste0('(',round(totalHeight,3)*100,'.0%)')),cex.lab=1.5,ylim=c(-0.25,2.25),xlim=c(-0.25,2.25),cex=2)
points(df$metline[which(df$stype=='Paedomorphwinter')],df$metheight[which(df$stype=='Paedomorphwinter')],pch=16,col='royalblue1')
lines(c(-0.5,3),c(0,0),lty='dashed',col='lightgrey')
lines(c(0,0),c(-1,3),lty='dashed',col='lightgrey')
lines(c(-10*ff[1],10*ff[1]),c(-10*ff[2],10*ff[2]))
lines(c(20*Mx,100*Lx),c(20*My,100*Ly))
points(mean(df$metline[which(df$stype=='Larvae')]),mean(df$metheight[which(df$stype=='Larvae')]),pch=16,col='seagreen1',cex=2)
points(mean(df$metline[which(df$stype=='Metamorph')]),mean(df$metheight[which(df$stype=='Metamorph')]),pch=16,col='coral4',cex=2)
points(mean(df$metline[which(df$stype=='Larvae')]),mean(df$metheight[which(df$stype=='Larvae')]),cex=2)
points(mean(df$metline[which(df$stype=='Metamorph')]),mean(df$metheight[which(df$stype=='Metamorph')]),cex=2)
points(mean(df$metline[which(df$stype=='Paedomorphwinter')]),mean(df$metheight[which(df$stype=='Paedomorphwinter')]),pch=16,col='royalblue1',cex=2)
points(mean(df$metline[which(df$stype=='Paedomorphwinter')]),mean(df$metheight[which(df$stype=='Paedomorphwinter')]),cex=2)





