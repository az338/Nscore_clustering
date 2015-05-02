

library(knitr)
library(xtable)

library(plyr)
library(reshape2)
## library(RColorBrewer)
## library(gridExtra)
library(gplots)
library(ggplot2)
library(GGally)
library(WriteXLS)
library(stringr)

# set working directory
#setwd('C:/Users//Azedine/Dropbox/ucc_az/bioseek_Xitong_Ellen/') # Windows path
#setwd('/home/az338/bioseek_Xitong_Ellen/NScore_cliques/') #Calculon path
setwd('/home/az338/Dropbox/ucc_az/bioseek_Xitong_Ellen/NScore_cliques/') # Rover path

# clusters we will focus the analysis on
clusters = c(18,20,22,25,37,45,72,17,39,60,62)



## Produce heatmap based on NScore
## cliques and averaged rank

# load target prediction ranks (incl. inactives)
load('NSPred.Rdata')

# load chem dataset
chems=read.table('smiles_NScore_standardized.tab',h=T,sep='\t',comment.char='',quote='"')

# get smiles which could not be loaded
error_smiles=system('more error_struct.txt | grep structure',intern=T)
error_smiles = unlist(llply(strsplit(error_smiles,'Cannot read structure : '),function(x) x[[2]]))
error_smiles = str_trim(error_smiles)

# TO BE PRINTED IN THE REPORT (+ mention 6 mol would not load)
#print(xtable(chems[match(error_smiles,chems$Structure),-c(1,4,6,7)],caption='Chemicals which cannot be processed by the target prediction tool and where hence removed from the analysis.',label='rmChems'),include.rownames=F)

# chemicals in the analysis
chems=chems[-match(error_smiles,chems$Structure),]


# generate proper data sturctures for ggplot
colnames(NScore_pred) = c('TARGET','UNIPROT',as.character(chems$SPID),'AVG.RANK')


# target pred dataset
targetPred = melt(NScore_pred[,-ncol(NScore_pred)],id=c('TARGET','UNIPROT'))
colnames(targetPred) = c('TARGET','UNIPROT','COMPOUND','RANK')
targetPred$CLUSTER = chems[match(targetPred$COMPOUND,chems$SPID),'CLUSTER_ID']


targetPred = subset(targetPred,CLUSTER  %in% clusters)

# Function to display the top 50 targets s
# for each chemicals in the clique
## d_ply(targetPred, .(CLUSTER), function(df) {
##     rankMat = acast(df,TARGET~COMPOUND,value.var='RANK')
##     l_ply(1:ncol(rankMat), function(j) {
##         head(rankMat[order(rankMat[,j]),j],50)
##     },.print=T)
## },.print=T)

d_ply(targetPred, .(CLUSTER), function(df) {
    file = file.path('/home/az338/Dropbox/ucc_az/bioseek_Xitong_Ellen/NScore_cliques',paste('Clique_',unique(df$CLUSTER),'_top50_targets.xlsx',sep=''))
    l = dlply(df, .(COMPOUND), function(df2) {
        return(head(df2[order(df2$RANK),c('TARGET','UNIPROT','RANK')],50))
    })
    cat(unique(as.character(df$COMPOUND)),'\n')
    WriteXLS('l',file,SheetNames=unique(as.character(df$COMPOUND)))
})

