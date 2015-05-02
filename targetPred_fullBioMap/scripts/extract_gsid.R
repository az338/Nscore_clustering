library(plyr)
library(RMySQL)
library(gdata)
library(xlsx)

# run on calculon

# set My SQL connection
drv = dbDriver("MySQL")
pass = readLines("~/.mysqlcred")
con = dbConnect(drv,user='az338',pass=pass,dbname='toxcast')

# load full BIOMAP compound list
cmp = read.csv('~/bioseek_Xitong_Ellen/BSK_EPA_AGENTS.csv',stringsAsFactor=F)


# for each compound in the file
# extract corresponding GSID from
# toxcast db
ids = llply(cmp$SPID, function(spid) {
    return(dbGetQuery(con,paste("select sa_sample_id,sa_gsid from sample where sa_sample_id='",spid,"'",sep='')))
})

# store chemical not in toxcast db
cmp_nodb=cmp[which(unlist(llply(ids,length)) == 0),]

# remove from initial dataset
cmp = cmp[-which(unlist(llply(ids,length)) == 0),]

# combine other results in data frame
ids = do.call('rbind',ids)
ids = subset(ids,!duplicated(ids))



# load ToxCast structures file
#tox21_struct=read.xls("/home/az338/data/toxcast_2014/ChemicalFiles/TOX21S_v4b_8599_23Oct2014.xlsx",quote='')
#tox21_struct=read.xlsx("/home/az338/data/toxcast_2014/ChemicalFiles/TOX21S_v4b_8599_23Oct2014.xlsx",sheetName = 'TOX21S_v4b_20141023')
tox21_struct=read.csv("/home/az338/data/toxcast_2014/ChemicalFiles/TOX21S_STRUCTURES.csv")
smiles = as.character(tox21_struct[match(ids$sa_gsid,tox21_struct$DSSTox_GSID),'STRUCTURE_SMILES']) 

# save chemicals with no structure in the structure file in ToxCast
write.table(rbind(cmp[which(is.na(smiles)),],cmp_nodb),"~/bioseek_Xitong_Ellen/NO_SMILES_fullBIOMAP.txt",row.names=F,col.names=T,quote=F,sep='\t')

# remove these chemicals from full dataset
idx = which(is.na(smiles))
full.dtset = data.frame(cmp,ids$sa_gsid,smiles)
full.dtset = full.dtset[-idx,]

# export smiles
write.table(full.dtset,'~/bioseek_Xitong_Ellen/smiles_fullBIOMAP.tab',sep='\t',col.names=T,row.names=F,quote=F)

# export smiles only file
writeLines(as.character(full.dtset$smiles),'~/bioseek_Xitong_Ellen/smiles_only_fullBIOMAP.smi')


