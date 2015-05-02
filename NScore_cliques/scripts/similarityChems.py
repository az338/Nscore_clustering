# libraries
from rdkit import Chem, DataStructs
import cPickle, gzip, sys, os, os.path
from collections import defaultdict
from optparse import OptionParser 
from rdkit.Chem import AllChem

# read chemicals
chem_f = open('/home/az338/bioseek_Xitong_Ellen/NScore_cliques/smiles_NScore_standardized.tab','r')
chems = []
header = chem_f.readline().split('\t')
for c in chem_f :
    chems.append(c.split('\t'))

chem_f.close()

# compute morgan Fingerprints radius 2 (2048bit)
header.append('Fingerprint')
allFP = []
for c in chems :
    try :
        mol = Chem.MolFromSmiles(c[0].strip('"'))
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        c.append(fp)
        allFP.append(fp)
    except :
        c.append(None)
        print 'Error for chemical '+ c[1] + '\nSkipping...'
        continue

# compute Tanimoto similarities 
Tanimoto_Simil = [] 
for c in chems :
    if(c[-1] != None) :
        Tanimoto_Simil.append(DataStructs.BulkCosineSimilarity(c[-1],allFP))
        
# write similarities to file
f = open('/home/az338/bioseek_Xitong_Ellen/NScore_cliques/tanimoto_scores.txt','w')
for s_vect,c in zip(Tanimoto_Simil,chems) :
    s_line = c[1]
    for s in s_vect :
        s_line+='\t'+str(s)
    f.write(s_line+'\n')
f.close()
        


