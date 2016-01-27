library(plyr)
library(reshape2)
library(Rcpi)

DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az//bioseek_Xitong_Ellen/NScore_cliques//data/'

# load chem dataset
chems = read.table(file.path(DATA_DIR,'smiles_NScore.tab'),h=T,sep='\t',comment.char='',quote='"',stringsAsFactors=F)
duplicated_spid = chems[which(duplicated(chems$smiles)),'SPID']
duplicated_chemIdx = which(duplicated(chems$smiles))
chems = chems[which(!duplicated(chems$smiles)),]

# load BioMap dataset
biomap = read.csv(file.path(DATA_DIR,'BSK_TOXCAST_NSCORE_export.csv'),header=T)
biomap = subset(biomap,SPID %in% chems$SPID)  # filter to those that are in the biomap dataset

# compute chem similarity matrix
mol = readMolFromSmi(file.path(DATA_DIR,'smiles_NScore_ST_corrected.smi'),type='mol')
fps = extractDrugExtendedComplete(mol,depth=4)
fps = fps[-duplicated_chemIdx,]

# compute tanimoto similarity
sim = laply(1:nrow(fps), function(x) {
  return(laply(1:nrow(fps), function(y) {
    return(calcDrugFPSim(fps[x,],fps[y,],fptype='complete',metric='tanimoto'))    
  }))
})
rownames(sim) = colnames(sim) = chems$SPID
simdf = melt(sim)
colnames(simdf) = c('SPID1','SPID2','SIM')

# them compute activity median disimilarities per assay + tests
i = 1
pvals = ddply(biomap, .(ACID), function(df) { # loop over assays
  
  cat(i,'ACID:',unique(df$ACID),'\n')
  i <<- i + 1
  # compute median actity per compound in the assay
  # i.e average NSCORE over concentrations w/ medians
  med_act = dlply(df, .(SPID), function(df2) {
    return(median(df2$NSCORE)) 
  })
  # compute biosimilarity data frame
  actmelt = ldply(1:length(med_act), function(i) {
    return(ldply(1:length(med_act), function(j) {
      return(data.frame(SPID1=names(med_act)[i],SPID2 = names(med_act)[j], DIFF = abs(med_act[[i]] - med_act[[j]])))
    }))
  })
  # print(head(actmelt))
  #print(dim(actmelt))
  
  # sort data frame and remove duplicates
  l_ply(1:nrow(actmelt), function(i) {
    actmelt[i,1:2] <<- sort(actmelt[i,1:2])
  })
  actmelt = actmelt[which(actmelt$SPID1 != actmelt$SPID2),]
  actmelt = actmelt[!duplicated(actmelt),]  
  actmelt = actmelt[order(as.character(actmelt$SPID1),as.character(actmelt$SPID2)),]
  
  # merge with similarity data frame
  df = merge(simdf,actmelt,by=c('SPID1','SPID2'))
  
  # discretise chemical similarity 
  df$simGp = unlist(llply(df$SIM, function(sim) {
    if(sim < .3) 
      return('LOW')
    else if(sim >= .3 & sim < .75) 
      return('MED')
    else
      return('HIGH')
  }))
  # pvalue
  p1 = wilcox.test(subset(df,simGp == 'HIGH')$DIFF,subset(df,simGp == 'MED')$DIFF)$p.value
  p2 = wilcox.test(subset(df,simGp == 'LOW')$DIFF,subset(df,simGp == 'MED')$DIFF)$p.value
  diff1 = abs(median(subset(df,simGp == 'HIGH')$DIFF) - median(subset(df,simGp == 'MED')$DIFF))
  diff2 = abs(median(subset(df,simGp == 'LOW')$DIFF) - median(subset(df,simGp == 'MED')$DIFF))
  data.frame(dHM=diff1,dLM=diff2,pHM=p1,pLM=p2)
})
write.csv(pvals,file.path(DATA_DIR,'median_act_diffs_HLM_and_pvals.csv'),row.names=F)



