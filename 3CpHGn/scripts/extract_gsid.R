
library(plyr)
library(RMySQL)
library(gdata)

drv = dbDriver("MySQL")
pass = readLines("~/.mysqlcred")
con = dbConnect(drv,user='az338',pass=pass,dbname='toxcast')

## load cytotoxic BioSeek compounds
ctBS_chems = read.csv('/home/az338/bioseek_Xitong_Ellen/3CpHGn/3CpHGn.csv',h=T)

# for each compound in the file
# extract corresponding GSID
ids = llply(ctBS_chems$spid, function(spid) {
    return(dbGetQuery(con,paste("select sa_sample_id,sa_gsid from sample where sa_sample_id='",spid,"'",sep='')))
})
ids = do.call('rbind',ids)
ids = subset(ids,!duplicated(ids))

# load ToxCast structures file
tox21_struct=read.xls("/home/az338/data/toxcast_2014/ChemicalFiles/TOX21S_v4b_8599_23Oct2014.xlsx")
smiles = as.character(tox21_struct[match(ids$sa_gsid,tox21_struct$DSSTox_GSID),'STRUCTURE_SMILES'])
compounds = tox21_struct[match(ids$sa_gsid,tox21_struct$DSSTox_GSID),c('DSSTox_GSID','TS_ChemName','TS_Description','STRUCTURE_SMILES')]

# export smiles
write.table(data.frame(ids,smiles),'~/bioseek_Xitong_Ellen/3CpHGn/smiles_3CpHGn.tab',sep='\t',col.names=T,row.names=F,quote=F)

# export chemInfo
save(compounds,file='~/bioseek_Xitong_Ellen/chemicals.Rdata')
