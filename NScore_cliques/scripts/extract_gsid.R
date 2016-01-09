
library(plyr)
library(RMySQL)
library(gdata)

drv = dbDriver("MySQL")
pass = readLines("~/.mysqlcred")
con = dbConnect(drv,user='az338',pass=pass,dbname='toxcast')

## load NScore BioSeek compounds
NScore_chems = read.csv('/scratch/az338/ucc-fileserver/ucc_az/bioseek_Xitong_Ellen/NScore_cliques/data/CliqueShedClusterAnnotation.csv',h=T)

# for each compound in the file
# extract corresponding GSID
ids = llply(NScore_chems$SPID, function(spid) {
    return(dbGetQuery(con,paste("select sa_sample_id,sa_gsid from sample where sa_sample_id='",spid,"'",sep='')))
})
ids = do.call('rbind',ids)
ids = subset(ids,!duplicated(ids))

# load ToxCast structures file
tox21_struct=read.xls("/scratch/az338/ucc-fileserver/ucc_az/toxCast_lincs_integration/data/toxcast_Oct14/ChemicalFiles/TOX21S_v4b_8599_23Oct2014.xlsx")
smiles = as.character(tox21_struct[match(ids$sa_gsid,tox21_struct$DSSTox_GSID),'STRUCTURE_SMILES']) 

# 3 compounds do not have a gsid in the tox21 structure file, and 1 (Clove leaf oil) is in the file but does not have a structure/smiles.
df = data.frame(NScore_chems,ids$sa_gsid,smiles)
df = df[-which(is.na(df$smiles)|df$smiles == ''),]

# also remove sodium chlorite as RDkit cannot process it
df = df[-grep('Sodium chlorite',df$AGENT_NAME),]

# export smiles
write.table(df,'/scratch/az338/ucc-fileserver/ucc_az/bioseek_Xitong_Ellen/NScore_cliques/data/smiles_NScore.tab',sep='\t',col.names=T,row.names=F,quote=F)
writeLines(df$smiles,'/scratch/az338/ucc-fileserver/ucc_az/bioseek_Xitong_Ellen/NScore_cliques/data/smiles_NScore.smi')


