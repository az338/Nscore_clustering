library(plyr)

# Read standardized smiles and corresponding
# cluster/clique number 

Chems = read.table('~/bioseek_Xitong_Ellen/NScore_cliques/smiles_NScore_standardized.tab',sep='\t',h=T)

# Run target prediction tool clique-wise (incl. negatives)
cliquePred = dlply(Chems,.(CLUSTER_ID), function(chems) {
    cat('Running PIDGIN for cluster num.',unique(chems$CLUSTER_ID),'with',nrow(chems),'chemicals\n')
    # write smiles in a file
    writeLines(as.character(chems$Structure),file.path('~/bioseek_Xitong_Ellen/NScore_cliques/cliqueTargets.csv'))
    # set working dir to PIDGIN dir
    setwd('~/PIDGIN')
    # run PIDGIN
    system('python predict_binary_heat.py a ~/bioseek_Xitong_Ellen/NScore_cliques/cliqueTargets.csv')
    # read and return rank matrix
    return(try(read.table('out_results_binary_heat.txt',h=F,sep="\t",comment.char="",quote="")))

    #remove output files
    system('rm ~/bioseek_Xitong_Ellen/NScore_cliques/cliqueTargets.csv')
    system('rm out_results_binary_heat.txt')
})
save(cliquePred, file = '~/bioseek_Xitong_Ellen/NScore_cliques/cliquePred.Rdata')



# Run target prediction tool clique-wise (single model)
cliquePredSingle = dlply(Chems,.(CLUSTER_ID), function(chems) {
    cat('Running PIDGIN for cluster num.',unique(chems$CLUSTER_ID),'with',nrow(chems),'chemicals\n')
    # write smiles in a file
    writeLines(as.character(chems$Structure),file.path('~/bioseek_Xitong_Ellen/NScore_cliques/cliqueTargets.csv'))
    # set working dir to PIDGIN dir
    setwd('~/PIDGIN/singlemodel')
    # run PIDGIN
    system('python predict_singlemodel.py a ~/bioseek_Xitong_Ellen/NScore_cliques/cliqueTargets.csv')
    # read and return rank matrix
    return(try(read.table('out_results_singlemodel.txt',h=F,sep="\t",comment.char="",quote="")))

    #remove output files
    system('rm ~/bioseek_Xitong_Ellen/NScore_cliques/cliqueTargets.csv')
    system('rm out_results_singlemodel.txt')
})
save(cliquePredSingle, file = '~/bioseek_Xitong_Ellen/NScore_cliques/cliquePredSingle.Rdata')

