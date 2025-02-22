library('vegan')
library('phyloseq')
library('ape')
library('PERMANOVA')
library('pairwiseAdonis')
library('factoextra')

#Read biom file
data <- import_biom('ASV_Table.biom')

#Find OTU table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

metadata<-read.csv('Mole_MetadataAC4.csv',row.names=1, header= TRUE)

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]


temp_phylo<-phyloseq(OTU_biom,TAX)
Mphylo<-prune_samples(colnames(OTU_biom)[which(colnames(OTU_biom) %in% rownames(metadata))],temp_phylo)

M2phylo<-phyloseq(otu_table(Mphylo),tax_table(Mphylo),sample_data(metadata))

M2phylo = prune_taxa(assigned_taxa_seqs, M2phylo)

M2phylo = prune_taxa(taxa_names(M2phylo)[which(tax_table(M2phylo)[,1]=='k__Bacteria')],M2phylo)

set.seed(1)

rMphylo = rarefy_even_depth(M2phylo, rngseed=1, replace=F, verbose = TRUE)

colvec<-c("coral4", "royalblue1","lightsteelblue1","seagreen2")
type<-sample_data(rMphylo)$StypeS
types<-c('Metamorph','w_Paedomorph','s_Paedomorph','Larvae')
#Convert names of groups to numbers
no_type<-rep(1,length(type))
for (k in 1:length(types)){
  no_type[which(type==types[k])]<-k
}


euclidean<-vegdist(t(otu_table(rMphylo)),method='euclidean',upper=TRUE,diag=TRUE,binary=FALSE)
euclidean_matrix<-as.matrix(euclidean,labels=TRUE)
pcoa_euclidean<-pcoa(euclidean)$vectors[,1:2]

#Run PCA analysis
pca <- prcomp(t(otu_table(rMphylo)))

#Find variation explained by axes
variance <- (pca$sdev)^2
varPercent <- variance/sum(variance) * 100

#Find loadings on axes
loadings <- pca$rotation

top3PCA1loadings<-names(sort(abs(loadings[,1]),decreasing=TRUE)[1:5])
top3PCA2loadings<-names(sort(abs(loadings[,2]),decreasing=TRUE)[1:5])

print('Top PCA1 loadings')
namestopPCA1<-list()
for (k in 1:5){
  findme<-which(rownames(tax_table(rMphylo))==top3PCA1loadings[k])
  namestopPCA1<-c(namestopPCA1,list(tax_table(rMphylo)[findme,]))
}

print('Top PCA2 loadings')
namestopPCA2<-list()
for (k in 1:5){
  findme<-which(rownames(tax_table(rMphylo))==top3PCA2loadings[k])
  namestopPCA2<-c(namestopPCA2,list(tax_table(rMphylo)[findme,]))
  
}

ff<-get_pca_var(pca)
contribtopPCA1<-sort(ff$contrib[,1],decreasing=TRUE)[1:5]
contribtopPCA2<-sort(ff$contrib[,2],decreasing=TRUE)[1:5]


par(mar = c(5, 5,5, 5))
plot(pcoa_euclidean,type='n',cex.lab=1.4,xlab=paste('PC1',paste0('(',round(varPercent[1],2),'%)')),ylab=paste('PC2',paste0('(',round(varPercent[2],2),'%)')))
ordihull(pcoa_euclidean,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_euclidean, display = "sites", pch=16, col = colvec[no_type],cex=1)
points(pcoa_euclidean,  pch=1, col = "black",cex=1)

#Find the centroid (based on mean) of each treatment group in the NMDS plot
xmeans_euclidean<-c()    #Value of the centroid along the first NMDS axis
ymeans_euclidean<-c()    #Value of the centroid along the second NMDS axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_euclidean<-c(xmeans_euclidean,mean(pcoa_euclidean[pts,1]))
  ymeans_euclidean<-c(ymeans_euclidean,mean(pcoa_euclidean[pts,2]))
}

#Plot the centroid for each group overtop of your NMDS scatterplot
for (k in 1:length(types)){
  points(xmeans_euclidean[k],ymeans_euclidean[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_euclidean[k],ymeans_euclidean[k],col='black',pch=1,cex=2)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the euclidean, bray, horn, morisita-horn and jaccard indices
euclidean<-vegdist(t(otu_table(rMphylo)),method='euclidean',upper=TRUE,diag=TRUE, binary=FALSE)

#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
euclidean_matrix<-as.matrix(euclidean,labels=TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform Statistical Tests

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##########################################  PERMANOVA  ##############################################################################################################################

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova<-list(Data=t(otu_table(rMphylo)),D=euclidean_matrix,Coefficient="Other")

type<-sample_data(rMphylo)$StypeS

#Alternatively, group summer and winter paeds
#type[type=='w_Paedomorph']<-'Paedomorph'
#type[type=='s_Paedomorph']<-'Paedomorph'

#Perform PERMANOVA comparing treatment groups
permanova_by_group=PERMANOVA(inputpermanova, factor(type),nperm=1000)

#If PERMANOVA by site is significant, figure out which groups differ from which other groups
if (permanova_by_group$pvalue<0.05 && length(unique(type))>2){
  posthoc_permanova<-pairwise.adonis(euclidean,as.factor(type),p.adjust.m='BH',perm=1000)
}

