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

FourHmodified <- function(x, class_grouping, core_fraction, boot_no, sample_no, core_average_abundance = 0, core_max_abundance = 0, core_all_average_abundance = 0, core_all_max_abundance = 0, replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE,seed=NULL,dist='Jaccard') {
  
  rlang::englue("var: {{ core_fraction }}")
  rlang::englue("var: {{ boot_no }}")
  rlang::englue("var: {{ sample_no }}")
  rlang::englue("var: {{ replace_hosts }}")
  
  #Define the number of reads
  if (is.null(reads)==TRUE){
    reads<-min(sample_sums(x))
  }
  
  #Define the fraction of hosts a microbe must be present on to be considered 'core'
  corefactor<-core_fraction
  
  #Define vectors for the four models and sigma from each trial
  fraction_INTERSECTION<-c()
  fraction_UNION<-c()
  fraction_GAIN<-c()
  fraction_LOSS<-c()
  fraction_LOSS1<-c()
  fraction_LOSS2<-c()
  fraction_LOSS12<-c()
  fraction_sigma<-c()
  fraction_UNION_P1<-c()
  fraction_UNION_P2<-c()
  
  #Define the seed
  if (is.null(seed)==TRUE){
    seed<-sample(10000,1)
  }
  
  #If you are not going to rarefy at each bootstrap, rarefy at the beginning (this uses a single rarefaction for all bootstrap samples)
  if (rarefy_each_step==FALSE){
    rared<-phyloseq::rarefy_even_depth(x,verbose=FALSE,sample.size=reads)
  }
  
  #For each bootstrap sample...
  for (m in 1:boot_no){
    
    if (rarefy_each_step==TRUE){
      #Rarefy the OTU table to the lowest read depth across all samples, with the option of specifying a lower read depth if you want to standardize across studies
      rared<-phyloseq::rarefy_even_depth(x,verbose=FALSE,sample.size=reads)
    }
    
    #Convert the OTU table to a matrix
    OTU1 = as(phyloseq::otu_table(rared), "matrix")
    
    #Create matrices for each host class (progenitor 1 == 1, hybrid == 2, progenitor 2 == 3)
    OTU1P1<-OTU1[,which(class_grouping==1)]
    OTU1H<-OTU1[,which(class_grouping==2)]
    OTU1P2<-OTU1[,which(class_grouping==3)]
    
    #Subsample from the matrices for each host class
    P1sub<-sample(1:dim(OTU1P1)[2],sample_no,replace=replace_hosts)
    Hsub<-sample(1:dim(OTU1H)[2],sample_no,replace=replace_hosts)
    P2sub<-sample(1:dim(OTU1P2)[2],sample_no,replace=replace_hosts)
    
    #Subsample matrices for this particular bootstrap
    OTU1P1sub<-OTU1P1[,P1sub]
    OTU1Hsub<-OTU1H[,Hsub]
    OTU1P2sub<-OTU1P2[,P2sub]
    OTU1sub<-cbind(OTU1P1,OTU1H,OTU1P2)
    
    #Find the core microbiome for each host class
    coreP1<-intersect(intersect(intersect(intersect(which(rowSums(sign(OTU1P1sub))>=round(corefactor*dim(OTU1P1sub)[2])), which(rowSums(OTU1P1sub)/(dim(OTU1P1sub)[2]*colSums(OTU1P1sub)[1])>=core_average_abundance)), which(apply(OTU1P1sub,1,max)/colSums(OTU1P1sub)[1]>=core_max_abundance)), which(rowSums(OTU1sub)/(dim(OTU1sub)[2]*colSums(OTU1sub)[1])>=core_all_average_abundance)), which(apply(OTU1sub,1,max)/colSums(OTU1sub)[1]>=core_all_max_abundance))
    coreH<-intersect(intersect(intersect(intersect(which(rowSums(sign(OTU1Hsub))>=round(corefactor*dim(OTU1Hsub)[2])), which(rowSums(OTU1Hsub)/(dim(OTU1Hsub)[2]*colSums(OTU1Hsub)[1])>=core_average_abundance)), which(apply(OTU1Hsub,1,max)/colSums(OTU1Hsub)[1]>=core_max_abundance)), which(rowSums(OTU1sub)/(dim(OTU1sub)[2]*colSums(OTU1sub)[1])>=core_all_average_abundance)), which(apply(OTU1sub,1,max)/colSums(OTU1sub)[1]>=core_all_max_abundance))
    coreP2<-intersect(intersect(intersect(intersect(which(rowSums(sign(OTU1P2sub))>=round(corefactor*dim(OTU1P2sub)[2])), which(rowSums(OTU1P2sub)/(dim(OTU1P2sub)[2]*colSums(OTU1P2sub)[1])>=core_average_abundance)), which(apply(OTU1P2sub,1,max)/colSums(OTU1P2sub)[1]>=core_max_abundance)), which(rowSums(OTU1sub)/(dim(OTU1sub)[2]*colSums(OTU1sub)[1])>=core_all_average_abundance)), which(apply(OTU1sub,1,max)/colSums(OTU1sub)[1]>=core_all_max_abundance))
    
    #Find the names of the core microbial taxa for each host class
    namesP1<-rownames(OTU1P1sub)
    namesH<-rownames(OTU1Hsub)
    namesP2<-rownames(OTU1P2sub)
    
    corenamesP1<-namesP1[coreP1]
    corenamesH<-namesH[coreH]
    corenamesP2<-namesP2[coreP2]
    
    #If you are using the Jaccard inspired 4H Index...
    if (dist == 'Jaccard'){
      #Intersection is the intersection between P1, P2 and H, divided by the total number of unique core microbial taxa across all three host classes
      fraction_INTERSECTION<-c(fraction_INTERSECTION,length(intersect(intersect(corenamesP1,corenamesP2),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #Union is the intersection between H and the microbes on P1 and P2 but not both, divided by the total number of unique core microbial taxa. To find the microbes on P1 and P2 but not both, we find the total list of unique microbes on P1 and P2 and remove those microbes shared by P1 and P2
      fraction_UNION<-c(fraction_UNION,(length(intersect(setdiff(unique(c(corenamesP1,corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #Gain is microbes on H that are not found on either P1 or P2, divided by the total number of unique core microbial taxa. To find the microbes found on H but not P1 or P2, we find the total number of microbes on H, and subtract those that intersect with the list of unique microbes found on P1 and P2
      fraction_GAIN<-c(fraction_GAIN,(length(corenamesH)-length(intersect(unique(c(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #Loss is the microbes found on P1, P2 or both that are not found on H, divided by the total number of unique core microbial taxa across all three host classes.
      fraction_LOSS<-c(fraction_LOSS,length(setdiff(unique(c(corenamesP1,corenamesP2)),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      fraction_LOSS12<-c(fraction_LOSS12,length(setdiff(intersect(corenamesP1,corenamesP2),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      fraction_LOSS1<-c(fraction_LOSS1,length(setdiff(setdiff(corenamesP1,corenamesP2),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      fraction_LOSS2<-c(fraction_LOSS2,length(setdiff(setdiff(corenamesP2,corenamesP1),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #sigma is the intersection of P1 and P2 divided by the total number of unique microbes found on P1 and P2
      fraction_sigma<-c(fraction_sigma,length(intersect(corenamesP1,corenamesP2))/length(unique(c(corenamesP1,corenamesP2))))
      fraction_UNION_P1<-c(fraction_UNION_P1,(length(intersect(setdiff(unique(c(corenamesP1)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      fraction_UNION_P2<-c(fraction_UNION_P2,(length(intersect(setdiff(unique(c(corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
    }
    else if (dist == 'Sorensen'){
      #Intersection is 3x the intersection between P1, P2 and H, divided by the total number of core microbial taxa on each host, summed across all three host classes
      fraction_INTERSECTION<-c(fraction_INTERSECTION,3*length(intersect(intersect(corenamesP1,corenamesP2),corenamesH))/length((c(corenamesH,corenamesP1,corenamesP2))))
      #Union is 2x the intersection between H and the microbes on P1 and P2 but not both, divided by the total number core microbial taxa on each host, summed across all three host classes. To find the microbes on P1 and P2 but not both, we find the total list of unique microbes on P1 and P2 and remove those microbes shared by P1 and P2
      fraction_UNION<-c(fraction_UNION,(2*length(intersect(setdiff(unique(c(corenamesP1,corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
      #Gain is microbes on H that are not found on either P1 or P2, divided by the total number core microbial taxa on each host, summed across all three host classes. To find the microbes found on H but not P1 or P2, we find the total number of microbes on H, and subtract those that intersect with the list of unique microbes found on P1 and P2
      fraction_GAIN<-c(fraction_GAIN,(length(corenamesH)-length(intersect(unique(c(corenamesP1,corenamesP2)),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
      #Loss is 1x the microbes found on P1 and/or P2 that are not found on H plus another 1x the microbes found on both P1 and P2, divided by the total number core microbial taxa on each host, summed across all three host classes.
      fraction_LOSS<-c(fraction_LOSS,(length(setdiff(unique(c(corenamesP1,corenamesP2)),corenamesH))+length(setdiff(intersect(corenamesP1,corenamesP2),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
      fraction_LOSS12<-c(fraction_LOSS12,2*length(setdiff(intersect(corenamesP1,corenamesP2),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      fraction_LOSS1<-c(fraction_LOSS1,length(setdiff(setdiff(corenamesP1,corenamesP2),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      fraction_LOSS2<-c(fraction_LOSS2,length(setdiff(setdiff(corenamesP2,corenamesP1),corenamesH))/length(unique(c(corenamesH,corenamesP1,corenamesP2))))
      #sigma is 2x the intersection of P1 and P2 divided by the total number of microbes found on P1 and P2 summed across P1 and P2
      fraction_sigma<-c(fraction_sigma,2*length(intersect(corenamesP1,corenamesP2))/length((c(corenamesP1,corenamesP2))))
      fraction_UNION_P1<-c(fraction_UNION_P1,(2*length(intersect(setdiff(unique(c(corenamesP1)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
      fraction_UNION_P2<-c(fraction_UNION_P2,(2*length(intersect(setdiff(unique(c(corenamesP2)),intersect(corenamesP1,corenamesP2)),corenamesH)))/length((c(corenamesH,corenamesP1,corenamesP2))))
    }
    else{
      print('Unrecognized model inspiration: Please use Jaccard or Sorensen')
    }
    
  }
  
  
  #Make a data frame of fractions
  fractiondf<-data.frame(fraction_INTERSECTION,fraction_GAIN,fraction_LOSS,fraction_LOSS1,fraction_LOSS2,fraction_LOSS12,fraction_UNION,fraction_UNION_P1,fraction_UNION_P2)
  
  
  return(fractiondf)
}


FourHmodifiedA <- function(x, class_grouping, core_fraction, boot_no, sample_no, core_average_abundance = 0, core_max_abundance = 0, core_all_average_abundance = 0, core_all_max_abundance = 0, replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE,seed=NULL,dist = 'Bray-Curtis',representative = 'mean',rescale_core = FALSE, use_microViz = 'no') {
  
  rlang::englue("var: {{ core_fraction }}")
  rlang::englue("var: {{ boot_no }}")
  rlang::englue("var: {{ sample_no }}")
  rlang::englue("var: {{ replace_hosts }}")
  
  
  if (use_microViz == 'yes'){
    #Check to make sure there is sample data and if not create a sample data dataframe with sample names
    x<-microViz::phyloseq_validate(x,verbose=FALSE)
  }
  
  
  #Define the number of reads
  if (is.null(reads)==TRUE){
    reads<-min(sample_sums(x))
  }
  
  #Define the fraction of hosts a microbe must be present on to be considered 'core'
  corefactor<-core_fraction
  
  #Define vectors for the four models and sigma from each trial
  fraction_INTERSECTION<-c()
  fraction_UNION<-c()
  fraction_GAIN<-c()
  fraction_LOSS<-c()
  fraction_LOSS1<-c()
  fraction_LOSS2<-c()
  fraction_LOSS12<-c()
  fraction_sigma<-c()
  fraction_UNION_P1<-c()
  fraction_UNION_P2<-c()
  
  #Define the seed
  if (is.null(seed)==TRUE){
    seed<-sample(10000,1)
  }
  
  #Put the class grouping into the phyloseq object
  sample_data(x)$classes<-class_grouping
  
  #If you are not going to rarefy at each bootstrap, rarefy at the beginning (this uses a single rarefaction for all bootstrap samples)
  if (rarefy_each_step==FALSE){
    rared<-phyloseq::rarefy_even_depth(x,verbose=FALSE,sample.size=reads)
  }
  
  
  #For each bootstrap sample...
  for (m in 1:boot_no){
    
    if (rarefy_each_step==TRUE){
      #Rarefy the OTU table to the lowest read depth across all samples, with the option of specifying a lower read depth if you want to standardize across studies
      rared<-phyloseq::rarefy_even_depth(x,verbose=FALSE,sample.size=reads)
    }
    
    
    #Make a phyloseq object for each host class
    rared1<-phyloseq::prune_samples(sample_data(rared)$classes == 1,rared)
    rared2<-phyloseq::prune_samples(sample_data(rared)$classes == 2,rared)
    rared3<-phyloseq::prune_samples(sample_data(rared)$classes == 3,rared)
    
    #Pick the individuals from each host class for this particular bootstrap
    pick1<-phyloseq::sample_names(rared1)[sample(1:length(sample_names(rared1)),sample_no,replace=replace_hosts)]
    pick2<-phyloseq::sample_names(rared2)[sample(1:length(sample_names(rared2)),sample_no,replace=replace_hosts)]
    pick3<-phyloseq::sample_names(rared3)[sample(1:length(sample_names(rared3)),sample_no,replace=replace_hosts)]
    
    #Make phyloseq objects for each host class for this particular bootstrap
    otu1<-phyloseq::prune_samples(pick1,rared1)
    otu2<-phyloseq::prune_samples(pick2,rared2)
    otu3<-phyloseq::prune_samples(pick3,rared3)
    
    #Find the names of the core microbial taxa for each host class for this particular bootstrap
    core1<-rownames(phyloseq::otu_table(otu1))[which(rowSums(sign(phyloseq::otu_table(otu1)))>=round(corefactor*sample_no))]
    core2<-rownames(phyloseq::otu_table(otu2))[which(rowSums(sign(phyloseq::otu_table(otu2)))>=round(corefactor*sample_no))]
    core3<-rownames(phyloseq::otu_table(otu3))[which(rowSums(sign(phyloseq::otu_table(otu3)))>=round(corefactor*sample_no))]
    
    #Provided a host class has at least one core microbial taxon, calculate the average or median number of reads associated with each core taxon for each host class
    #Progenitor 1
    if (length(core1)>0){
      coreotu1<-phyloseq::prune_taxa(core1,otu1)
      if (representative == 'mean'){
        mean1<-rowMeans(phyloseq::otu_table(coreotu1))
      }
      else if (representative == 'median'){
        mean1<-apply(phyloseq::otu_table(coreotu1),1,median,na.rm=TRUE)/sum(apply(phyloseq::otu_table(coreotu1),1,median,na.rm=TRUE))
      }
      else{
        print('Method for finding representative microbial abundances not recognized: Please use either mean or median')
      }
    }
    #Hybrid
    if (length(core2)>0){
      coreotu2<-phyloseq::prune_taxa(core2,otu2)
      if (representative == 'mean'){
        mean2<-rowMeans(phyloseq::otu_table(coreotu2))
      }
      else if (representative == 'median'){
        mean1<-apply(phyloseq::otu_table(coreotu2),1,median,na.rm=TRUE)/sum(apply(phyloseq::otu_table(coreotu2),1,median,na.rm=TRUE))
      }
      else{
        print('Method for finding representative microbial abundances not recognized: Please use either mean or median')
      }
    }
    #Progenitor 2
    if (length(core3)>0){
      coreotu3<-phyloseq::prune_taxa(core3,otu3)
      if (representative == 'mean'){
        mean3<-rowMeans(phyloseq::otu_table(coreotu3))
      }
      else if (representative == 'median'){
        mean1<-apply(phyloseq::otu_table(coreotu3),1,median,na.rm=TRUE)/sum(apply(phyloseq::otu_table(coreotu3),1,median,na.rm=TRUE))
      }
      else{
        print('Method for finding representative microbial abundances not recognized: Please use either mean or median')
      }
    }
    
    #If you are going to rescale the core microbiome reads to total one, divide the average number of reads for each taxon by the total average number of reads across all taxa in the core
    if (rescale_core == TRUE){
      if (length(core1)>0){
        mean1<-mean1/sum(mean1)
      }
      if (length(core2)>0){
        mean2<-mean2/sum(mean2)
      }
      if (length(core3)>0){
        mean3<-mean3/sum(mean3)
      }
      
    }
    
    #If there is no core microbiome on any host class, set all parameters to zero
    if (length(core1) == 0 && length(core2) == 0 && length(core3) == 0){
      print('no core microbiome on any class')
      n_intersection<-0
      n_union<-0
      n_union1<-0
      n_union2<-0
      n_gain<-0
      n_loss<-0
    }
    
    #As long as at least one core microbial taxon on at least one host class...
    else{
      
      #If progenitor 1 has no core microbial taxa, arbitrarly pick a host class that DOES have a core and use those microbial taxon names... but set their abundances to zero
      if (length(core1)==0){
        if (length(core2)>0){
          mean1<-mean2
          mean1[1:length(mean1)]<-0
        }
        else{
          mean1<-mean3
          mean1[1:length(mean1)]<-0
        }
      }
      
      #If the hybrid has no core microbial taxa, arbitrarly pick a host class that DOES have a core and use those microbial taxon names... but set their abundances to zero
      if (length(core2)==0){
        if (length(core1)>0){
          mean2<-mean1
          mean2[1:length(mean2)]<-0
        }
        else{
          mean2<-mean3
          mean2[1:length(mean2)]<-0
        }
      }
      
      #If progenitor 3 has no core microbial taxa, arbitrarly pick a host class that DOES have a core and use those microbial taxon names... but set their abundances to zero
      if (length(core3)==0){
        if (length(core2)>0){
          mean3<-mean2
          mean3[1:length(mean3)]<-0
        }
        else{
          mean3<-mean1
          mean3[1:length(mean3)]<-0
        }
      }
      
      
      #make a matrix with the core microbial taxa and abundances from progenitor 1 and the hybrid
      vv<-merge(mean1,mean2,by='row.names',all=TRUE)
      vv[which(is.na(vv[,2])),2]<-0
      vv[which(is.na(vv[,3])),3]<-0
      rownames(vv)<-vv[,1]
      vv<-vv[,2:3]
      #make a matrix with the microbial taxa and abundances from progenitor 1, the hybrid and progenitor 2
      ww<-merge(vv,mean3,by='row.names',all=TRUE)
      ww[which(is.na(ww[,2])),2]<-0
      ww[which(is.na(ww[,3])),3]<-0
      ww[which(is.na(ww[,4])),4]<-0
      rownames(ww)<-ww[,1]
      ww<-ww[,2:4]
      colnames(ww)<-c('parent1','hybrid','parent2')
      
      #Calculate the number of reads on the progenitors to find sigma (this is regardless of the model that these reads fall into)
      #reads on both progenitors double counting the reads that they share
      sigma_total<-ww[,1]+ww[,3]
      #Only the reads that the two progenitors share
      sigma_shared<-apply(cbind(ww[,1],ww[,3]),1,FUN=min)
      #create a matrix for subtracting the shared progenitor 1/progenitor 2 reads from progenitor 1 and progenitor 2
      sigma_unique_matrix<-cbind(ww[,1],ww[,3])-sigma_shared
      #Only the reads that are on one or other of the progenitors, but not both
      sigma_unique<-sigma_unique_matrix[,1]+sigma_unique_matrix[,2]
      
      
      #for each microbial taxon, find the minimum number of reads across each host class... these are the reads shared by all three host classes for core taxa
      intersection<-apply(ww, 1, FUN = min)
      
      #subtract the shared reads from the matrix
      ww<-ww-intersection
      
      #For the reads that remain, consider progenitor 1 and the hybrid... for each microbial taxon, find the minimum number of reads across these two host classes... these are the reads shared by progenitor 1 and the hybrid for core taxa
      union1<-apply(ww[,1:2],1,FUN=min)
      #create a matrix for subtracting the progenitor 1/hybrid reads from progenitor 1 and the hybrid (subtract nothing from progenitor 2)
      union1m<-cbind(union1,union1,rep(0,length(union1)))
      ww<-ww-union1m
      
      #For the reads that remain, consider progenitor 2 and the hybrid... for each microbial taxon, find the minimum number of reads across these two host classes... these are the reads shared by progenitor 2 and the hybrid for core taxa
      union2<-apply(ww[,2:3],1,FUN=min)
      #create a matrix for subtracting the progenitor 2/hybrid reads from progenitor 2 and the hybrid (subtract nothing from progenitor 1)
      union2m<-cbind(rep(0,length(union2)),union2,union2)
      ww<-ww-union2m
      
      #The full union model are those reads shared with EITHER progenitor
      union<-union1+union2
      
      #Any reads remaining on the hybrid are attributed to the gain model
      gain<-ww[,2]
      
      #Any reads remaining on EITHER progenitor or both are attributed to the loss model
      #Loss double counting reads on both progenitors
      loss<-(ww[,1]+ww[,3])
      #Only the shared loss reads found on both progenitors
      loss_shared<-apply(cbind(ww[,1],ww[,3]),1,FUN=min)
      #create a matrix for subtracting the shared progenitor 1/progenitor 2 reads from progenitor 1 and progenitor 2
      loss_unique_matrix<-cbind(ww[,1],ww[,3])-loss_shared
      #Only the loss reads found on one or other progenitor, but not both
      loss_unique<-loss_unique_matrix[,1]+loss_unique_matrix[,2]
      
      
      
      #sum up the total reads for each model
      if (dist == 'Bray-Curtis'){
        #multiply the intersection reads by 3 (since there was one copy on EACH of the three host classes)
        s_intersection<-3*sum(intersection)
        #multiply the union reads by 2 (since there was one copy on a progenitor and one copy on the host)
        s_union<-2*sum(union)
        s_union1<-2*sum(union1)
        s_union2<-2*sum(union2)
        #do not multiply gain or loss by anything, since these were summed reads on individuals (or, if on both progenitors, were counted already because all remaining reads on progenitors were summed)
        s_gain<-sum(gain)
        s_loss<-sum(loss)
        s_loss1<-sum(loss_unique_matrix[,1])
        s_loss2<-sum(loss_unique_matrix[,2])
        s_loss12<-2*sum(loss_shared)
        #for sigma, calculate the number of shared reads relative to the number of total reads across parents
        s_sigma<-2*sum(sigma_shared)/sum(sigma_total)
      }
      else if (dist == 'Ruzicka'){
        #multiply the intersection reads by 3 (since there was one copy on EACH of the three host classes)
        s_intersection<-sum(intersection)
        #multiply the union reads by 2 (since there was one copy on a progenitor and one copy on the host)
        s_union<-sum(union)
        s_union1<-sum(union1)
        s_union2<-sum(union2)
        #do not multiply gain or loss by anything, since these were summed reads on individuals (or, if on both progenitors, were counted already because all remaining reads on progenitors were summed)
        s_gain<-sum(gain)
        s_loss<-sum(loss_shared+loss_unique)
        s_loss1<-sum(loss_unique_matrix[,1])
        s_loss2<-sum(loss_unique_matrix[,2])
        s_loss12<-sum(loss_shared)
        #for sigma, calculate the number of shared reads relative to the number of total reads across parents
        s_sigma<-sum(sigma_shared)/sum(sigma_shared+sigma_unique)
      }
      else{
        print('Unrecognized model inspiration: Please use Bray-Curtis or Ruzicka')
      }
      
      #turn the reads for each model into a fraction of the total reads
      n_intersection<-s_intersection/(s_intersection+s_union+s_gain+s_loss)
      n_union<-s_union/(s_intersection+s_union+s_gain+s_loss)
      n_union1<-s_union1/(s_intersection+s_union+s_gain+s_loss)
      n_union2<-s_union2/(s_intersection+s_union+s_gain+s_loss)
      n_gain<-s_gain/(s_intersection+s_union+s_gain+s_loss)
      n_loss<-s_loss/(s_intersection+s_union+s_gain+s_loss)
      n_loss1<-s_loss1/(s_intersection+s_union+s_gain+s_loss)
      n_loss2<-s_loss2/(s_intersection+s_union+s_gain+s_loss)
      n_loss12<-s_loss12/(s_intersection+s_union+s_gain+s_loss)
      
    }
    
    
    #Find the fraction of each 'model' represented by the hybrid microbiomes
    fraction_INTERSECTION<-c(fraction_INTERSECTION,n_intersection)
    fraction_UNION<-c(fraction_UNION,n_union)
    fraction_GAIN<-c(fraction_GAIN,n_gain)
    fraction_LOSS<-c(fraction_LOSS,n_loss)
    fraction_LOSS1<-c(fraction_LOSS1,n_loss1)
    fraction_LOSS2<-c(fraction_LOSS2,n_loss2)
    fraction_LOSS12<-c(fraction_LOSS12,n_loss12)
    fraction_sigma<-c(fraction_sigma,s_sigma)
    fraction_UNION_P1<-c(fraction_UNION_P1,n_union1)
    fraction_UNION_P2<-c(fraction_UNION_P2,n_union2)
    
    
  }
  
  
  #Make a data frame of fractions
  fractiondf<-data.frame(fraction_INTERSECTION,fraction_GAIN,fraction_LOSS,fraction_LOSS1,fraction_LOSS2,fraction_LOSS12,fraction_UNION,fraction_UNION_P1,fraction_UNION_P2)
  
  
  return(fractiondf)
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

M2phylo = prune_taxa(assigned_taxa_seqs, M2phylo)

M2phylo = prune_taxa(taxa_names(M2phylo)[which(tax_table(M2phylo)[,1]=='k__Bacteria')],M2phylo)

set.seed(1)

rMphylo = rarefy_even_depth(M2phylo, rngseed=1, replace=F, verbose = TRUE)

grouping<-sample_data(rMphylo)$StypeS
grouping[grouping=='Larvae']<-1
grouping[grouping=='s_Paedomorph']<-2
grouping[grouping=='w_Paedomorph']<-2
grouping[grouping=='Metamorph']<-3
grouping<-as.numeric(grouping)

#Calculate 4H-index Jaccard
ff<-FourHmodified(rMphylo,grouping,0.09,100,10,seed = 1)
mm<-round(100*colMeans(data.frame(ff)),2)

mmdf<-data.frame(c('MLP','P','ML','MP','L','LP','M'),c(mm[[1]],mm[[2]],mm[[6]],mm[[9]],mm[[4]],mm[[8]],mm[[5]]))
colnames(mmdf)<-c('group','FourH')

# Basic piechart
ggplot(mmdf, aes(x="", y=FourH, fill=group)) +
  geom_bar(stat="identity", width=1,fill=c('white','lightsteelblue1','lightsteelblue1','coral4','coral4','seagreen1',"seagreen1"),color='black') +
  coord_polar("y", start=0)+theme_void()

#Calculate 4H-index Bray-Curtis
gg<-FourHmodifiedA(rMphylo,grouping,0.09,100,10,seed=1)
nn<-round(100*colMeans(data.frame(gg)),2)

nndf<-data.frame(c('MLP','P','ML','MP','L','LP','M'),c(nn[[1]],nn[[2]],nn[[6]],nn[[9]],nn[[4]],nn[[8]],nn[[5]]))
colnames(nndf)<-c('group','FourH')

# Basic piechart
ggplot(nndf, aes(x="", y=FourH, fill=group)) +
  geom_bar(stat="identity", width=1,fill=c('white','lightsteelblue1','lightsteelblue1','coral4','coral4','seagreen1',"seagreen1"),color='black') +
  coord_polar("y", start=0)+theme_void()



