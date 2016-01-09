
library(plyr)


############################################################
## Run target prediction tool (incl. negatives)
############################################################
# write smiles in a file

# set working dir to PIDGIN dir
setwd('~/ucc-fileserver/ucc_az/PIDGIN')
# run PIDGIN for different threshold method
method = list('a','f','p','r')
system('rm ~/ucc-fileserver/ucc_az/bioseek_Xitong_Ellen/NScore_cliques/data/error_struct.txt')
l_ply(method, function(p) {
    system(paste('python predict_binary_heat.py ', p,' ~/ucc-fileserver/ucc_az/bioseek_Xitong_Ellen/NScore_cliques/data/smiles_NScore_ST_corrected.smi 2>> ~/ucc-fileserver/ucc_az/bioseek_Xitong_Ellen/NScore_cliques/data/error_struct.txt',sep=''))
# read and save rank matrix
    NScore_pred = read.table('out_results_binary_heat.txt',h=F,sep="\t",comment.char="",quote="",stringsAsFactors=F)
    save(NScore_pred, file = paste('~/ucc-fileserver/ucc_az/bioseek_Xitong_Ellen/NScore_cliques/data/NSPred_binary_heat_',p,'.Rdata',sep=''))
#remove output file
    system('rm out_results_binary_heat.txt')
})




#######################################     
## Run target prediction (single model)
#######################################
## # set working dir to PIDGIN dir
## setwd('~/PIDGIN/singlemodel')
## # run PIDGIN
## system('python predict_singlemodel_ranked_number.py ~/bioseek_Xitong_Ellen/NScore_cliques/SMILES_TMP.csv')
## # read and return rank matrix
## NScore_predSingle = read.table('out_results_singlemodel_ranked_number.txt',h=F,sep="\t",comment.char="",quote="")
## #remove input/output files
## system('rm ~/bioseek_Xitong_Ellen/NScore_cliques/SMILES_TMP.csv')
## system('rm out_results_singlemodel_ranked_number.txt')
## save(NScore_predSingle, file = '~/bioseek_Xitong_Ellen/NScore_cliques/NSPred_Single.Rdata')

