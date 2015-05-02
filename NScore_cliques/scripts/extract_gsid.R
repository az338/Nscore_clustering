
library(plyr)
library(RMySQL)
library(gdata)

drv = dbDriver("MySQL")
pass = readLines("~/.mysqlcred")
con = dbConnect(drv,user='az338',pass=pass,dbname='toxcast')

## load NScore BioSeek compounds
NScore_chems = read.csv('/home/az338/bioseek_Xitong_Ellen/NScore_cliques/CliqueShedClusterAnnotation.csv',h=T)

# for each compound in the file
# extract corresponding GSID
ids = llply(NScore_chems$SPID, function(spid) {
    return(dbGetQuery(con,paste("select sa_sample_id,sa_gsid from sample where sa_sample_id='",spid,"'",sep='')))
})
ids = do.call('rbind',ids)
ids = subset(ids,!duplicated(ids))

# load ToxCast structures file
tox21_struct=read.xls("/home/az338/data/toxcast_2014/ChemicalFiles/TOX21S_v4b_8599_23Oct2014.xlsx")
smiles = as.character(tox21_struct[match(ids$sa_gsid,tox21_struct$DSSTox_GSID),'STRUCTURE_SMILES']) 

# 3 compounds do not have a gsid in the tox21 structure file, and 1 (Clove leaf oil) is in the file but does not have a structure/smiles.

# export smiles
write.table(data.frame(NScore_chems,ids$sa_gsid,smiles),'~/bioseek_Xitong_Ellen/NScore_cliques/smiles_NScore.tab',sep='\t',col.names=T,row.names=F,quote=F)



