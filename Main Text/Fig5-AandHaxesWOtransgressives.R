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
library('StereoMorph')

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
    metamorphosis_position_P2<-(1-disttemp2[2,4]/disttemp2[1,2])*disttemp2[1,2]/disttemp[1,2]
  }else if (disttemp2[1,4]>disttemp2[2,4]){
    metamorphosis_position_P2<-(1+disttemp2[2,4]/disttemp2[1,2])*disttemp2[1,2]/disttemp[1,2]
  }else{
    metamorphosis_position_P2<-(0-disttemp2[1,4]/disttemp2[1,2])*disttemp2[1,2]/disttemp[1,2]
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
Pphylo<-prune_samples(sample_data(rMphylo)$StypeS %in% c('w_Paedomorph','s_Paedomorph'),rMphylo)
Pphylobig<-prune_samples(sample_data(Pphylo)$svl>51,Pphylo)
Pphylosmall<-prune_samples(sample_data(Pphylo)$svl<52,Pphylo)

Lcentroid<-as.vector(rowMeans(otu_table(Lphylo)))
Mcentroid<-as.vector(rowMeans(otu_table(Mphylo)))
Pcentroid<-as.vector(rowMeans(otu_table(Pphylobig)))
Ocentroid<-orthogonalProjectionToLine(Pcentroid, l1 = Lcentroid, l2 = Mcentroid)

locale<-rep(0,length(taxa_names(Lphylo)))
testme<-rep(0,length(taxa_names(Lphylo)))
dd<-data.frame(t(otu_table(rMphylo)[1,]),sample_data(rMphylo)$StypeS)
dd[,2][which(dd[,2]=='w_Paedomorph')]<-'Paedomorph'
dd[,2][which(dd[,2]=='s_Paedomorph')]<-'Paedomorph'

colnames(dd)<-c('count','type')
pv1<-c()
pv2<-c()
for (k in 1:length(locale)){
  dd[,1]<-t(otu_table(rMphylo)[k,])
  test<-pairwise.wilcox.test(dd$count, dd$type,p.adjust.method = "BH")
  pv1<-c(pv1,test$p.value[2,1])
  pv2<-c(pv2,test$p.value[2,2])
  if (is.nan(test$p.value[2,1])==FALSE && is.nan(test$p.value[2,2])==FALSE){
  if (max(test$p.value[2,])<0.05){
    locale[k]<-1
    testme[k]<-12
  } else if (test$p.value[2,1]<0.05){
    testme[k]<-1
  } else if (test$p.value[2,2]<0.05){
    testme[k]<-2
  }
  }else if (is.nan(test$p.value[2,1])==FALSE){
    if (test$p.value[2,1]<0.05){
      testme[k]<-1
    }
  }else if (is.nan(test$p.value[2,2])==FALSE){
    if (test$p.value[2,2]<0.05){
      testme[k]<-2
    }
  }
}


hitlist<-which(locale>0)

differentiate<-c()
for (k in 1:length(hitlist)){
  print(tax_table(rMphylo)[hitlist[k],5:6])
  Lmed<-mean(as.vector(otu_table(Lphylo)[hitlist[k],]))
  Pmed<-mean(as.vector(otu_table(Pphylo)[hitlist[k],]))
  Mmed<-mean(as.vector(otu_table(Mphylo)[hitlist[k],]))
  if (Pmed<Lmed && Pmed<Mmed){
    differentiate<-c(differentiate,-2)
  }else if (Pmed>Lmed && Pmed>Mmed){
    differentiate<-c(differentiate,2)
  }else if (Pmed>Lmed && Pmed<Mmed){
    differentiate<-c(differentiate,1)
  }else if (Pmed>Mmed && Pmed<Lmed){
    differentiate<-c(differentiate,1)
  }else{differentiate<-c(differentiate,0)}
}

ff<-taxa_names(rMphylo)[which(locale>0)]

rcutMphylo<-rMphylo
otu_table(rcutMphylo)[which(locale>0),]<-0

Mphylo<-prune_samples(sample_data(rcutMphylo)$StypeS=='Metamorph',rcutMphylo)
Lphylo<-prune_samples(sample_data(rcutMphylo)$StypeS=='Larvae',rcutMphylo)
Pphylo<-prune_samples(sample_data(rcutMphylo)$StypeS %in% c('w_Paedomorph','s_Paedomorph'),rcutMphylo)
Pphylobig<-prune_samples(sample_data(Pphylo)$svl>51,Pphylo)
Pphylosmall<-prune_samples(sample_data(Pphylo)$svl<52,Pphylo)

Aphylo<-merge_phyloseq(Lphylo,Mphylo,Pphylosmall)
wAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='w_Paedomorph',Aphylo)
sAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='s_Paedomorph',Aphylo)
lAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='Larvae',Aphylo)
mAphylo<-prune_samples(sample_data(Aphylo)$StypeS=='Metamorph',Aphylo)


Ametline<-c()
Ametheight<-c()
for (k in 1:length(sample_names(Aphylo))){
  temp<-Find_Projections(Lcentroid,Mcentroid,Pcentroid,Ocentroid,as.vector(otu_table(Aphylo)[,k]))
  Ametline<-c(Ametline,temp[1])
  Ametheight<-c(Ametheight,temp[2])
}



df<-data.frame(sample_data(Aphylo)$StypeS,Ametline,Ametheight)
colnames(df)<-c('stype','metline','metheight')

par(mar=c(5,5,5,2))
plot(df$metline,df$metheight,xlab='Ancestral (Metamorph) Axis',ylab='Heterokairic (Paedomorph) Axis',cex.lab=1.5,ylim=c(-0.1,1.9),xlim=c(-0.9,0.9))
points(df$metline[which(df$stype=='s_Paedomorph')],df$metheight[which(df$stype=='s_Paedomorph')],pch=16,col='lightsteelblue1')
points(df$metline[which(df$stype=='Metamorph')],df$metheight[which(df$stype=='Metamorph')],pch=16,col='coral4')
points(df$metline[which(df$stype=='Larvae')],df$metheight[which(df$stype=='Larvae')],pch=16,col='seagreen1')
points(df$metline[which(df$stype=='w_Paedomorph')],df$metheight[which(df$stype=='w_Paedomorph')],pch=16,col='royalblue1')
lines(c(-10,10),c(0,0),lty='dashed')
lines(c(0,0),c(-10,50),lty='dashed')

kw_maxis_4<-kruskal.test(metline ~ stype, data = df)
pw_maxis_4<-pairwise.wilcox.test(df$metline, df$stype,p.adjust.method = "BH")
kw_paxis_4<-kruskal.test(metheight ~ stype, data = df)
pw_paxis_4<-pairwise.wilcox.test(df$metheight, df$stype,p.adjust.method = "BH")

df$stype<-factor(df$stype,levels=c('Larvae','s_Paedomorph','w_Paedomorph','Metamorph'))


#Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
p_value_ns<-'no'
pshow<-'star'

##########################Metamorphosis Axis Violin Plot#########################################################################################################################################

#Define violin plot
p_maxis_4 <- ggplot(df, aes(x=stype, y=metline)) + geom_violin(aes(fill=stype)) + theme(axis.title.x = element_blank())+ylab('Ancestral (Metamorph) Axis')
#Choose the size of font for the axes titles and labels
p_maxis_4<-p_maxis_4+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 12))
#Choose the size of font for the legend title and lables
p_maxis_4<-p_maxis_4+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_maxis_4<-p_maxis_4+scale_fill_manual(values=c("seagreen2", "lightsteelblue1",'royalblue1' ,"coral4"))
#Add boxplots inside the violins
p_maxis_4<-p_maxis_4+geom_boxplot(aes(fill=stype),width=0.1)

#If there are significant differences in richness between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_maxis_4$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-1
  y_step<-0.16
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_maxis_4$p.value))){
    for (j in 1:k){
      if (rownames(pw_maxis_4$p.value)[k]!=colnames(pw_maxis_4$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_maxis_4$p.value[k,j]<0.1 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_maxis_4$p.value)[k])
          group2<-c(group2,colnames(pw_maxis_4$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_maxis_4$p.value[k,j])),5)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  #pstar<-pstar[c(2,6,4,3,5,1)]
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  stat.test_maxis_4<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_maxis_4<-p_maxis_4+stat_pvalue_manual(stat.test_maxis_4,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_maxis_4 #SI Fig. 1.1a


##########################Metamorphosis Axis Violin Plot#########################################################################################################################################

#Define violin plot
p_paxis_4 <- ggplot(df, aes(x=stype, y=metheight)) + geom_violin(aes(fill=stype)) + theme(axis.title.x = element_blank())+ylab('Heterokairic (Paedomorph) Axis')
#Choose the size of font for the axes titles and labels
p_paxis_4<-p_paxis_4+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 12))
#Choose the size of font for the legend title and lables
p_paxis_4<-p_paxis_4+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_paxis_4<-p_paxis_4+scale_fill_manual(values=c("seagreen2", "lightsteelblue1",'royalblue1' ,"coral4"))
#Add boxplots inside the violins
p_paxis_4<-p_paxis_4+geom_boxplot(aes(fill=stype),width=0.1)

#If there are significant differences in richness between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_paxis_4$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-0.12
  y_step<-0.02
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_paxis_4$p.value))){
    for (j in 1:k){
      if (rownames(pw_paxis_4$p.value)[k]!=colnames(pw_paxis_4$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_paxis_4$p.value[k,j]<0.1 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_paxis_4$p.value)[k])
          group2<-c(group2,colnames(pw_paxis_4$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_paxis_4$p.value[k,j])),5)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_paxis_4<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_paxis_4<-p_paxis_4+stat_pvalue_manual(stat.test_paxis_4,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_paxis_4 #SI Fig. 1.1a


intermediate<-length(which(differentiate==1))
lowside<-length(which(differentiate==-2))
highside<-length(which(differentiate==2))

dfs<-data.frame(c('intermediate','transgressive (below)','transgressive (above)'),c(intermediate,lowside,highside))
colnames(dfs)<-c('difference','amplicon sequence variants')

ggbarplot(dfs, "difference", "amplicon sequence variants",
          fill = 'difference', color = "difference",palette = c("purple", "blue", "red"),
          label = TRUE, lab.pos = "out", lab.col = "black",xlab = '',ylab='Amplicon Sequence Variants (ASVs)')+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 13))

