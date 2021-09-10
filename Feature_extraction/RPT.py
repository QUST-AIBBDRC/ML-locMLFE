import itertools
import random
import scipy.io as sio
import pickle as p
from itertools import combinations
import numpy as np
import math
#LG is some predefined constant
LG =  5
import math

def trigrams(matrix):
    return_matrix = []
    for x in range(len(matrix[0])):
      for y in range(len(matrix[0])):
        for z in range(len(matrix[0])):
          for i in range( len(matrix) - 2):
            value = matrix[i][x] * matrix[i+1][y] * matrix[i+2][z]
            return_matrix.append( value )
    return return_matrix

def residue_probing_transform(matrix):####RPT
  features = []
  #lets assume matrix[0] is 20 which is generally the case
  #sum of nth aa group in pssm
  matrix_sums = list(map(sum, zip(*matrix)))
  res_trnfr = []
  for i in range(len(matrix[0])):
    for j in range(len(matrix[0])):
      features.append( ( matrix_sums[i]+ matrix_sums[j]) / len(matrix) )
  return features

def evolutionary_diff(matrix, d = 30):####evolutionary distance transformation （EDT）
  features = []
  #assume all features ares same len and = 20
  for x in range(len(matrix[0] ) ):
    for y in range( len(matrix[0] ) ):
      value = 0
      for i in range( len(matrix) - d ):
        value += (matrix[i][x] - matrix[i+d][y]) ** 2
      value = value / (len(matrix) - d)
      features.append(value )
  #print(len(features) )
  return features
def pssm_sdt_func(matrix):
    #construct empty matrix of size 20 * LG
    pssm_sdt = [[0 for x in range( 20 )] for y in range(LG) ]
    #loop through all i's i.e amino acids
    for lg in range(0, LG):
        for i in range(len(matrix[0] ) ): #number of features
            pssm_sdt[lg][i] = 0
            for j in range(0, len(matrix) - lg):
                pssm_sdt[lg][i] += matrix[j][i] * matrix[j + lg][i] / (len(matrix) - lg)
    return pssm_sdt


#pssm distance transformation of different proteins
def pssm_ddt_func(matrix):
    #lets first calculate a permutations of all indexes in the matrix
    #we can calculate indices for pssm_ddt later on
    a = [x for x in range(20 )] # 20 is total number of amino acids
    #generates matrix of size 380
    indices_matrix = list(itertools.permutations(a, 2)) #combinations of two residues

    #empty matrix
    pssm_ddt = [[0 for x in range(19 * 20 )] for y in range(LG)]

    for lg in range(0, LG):
        for i1 in range(len(matrix[1])): #number of features
            for i2 in range(len(matrix[1])) :
                if i1 == i2:
                    continue
                index = indices_matrix.index((i1, i2))
                for j in range(0, len(matrix) - lg):
                    pssm_ddt[lg][index] += matrix[j][i1] * matrix[j+ lg][i2] / (len(matrix) - lg )
    return pssm_ddt

def feature_space(matrix, pssm_sdt=True, pssm_ddt=True):
    pssm_ddt1 = []
    pssm_sdt1 = []
    if pssm_ddt:
      pssm_ddt1 = [item for sublist in  pssm_ddt_func(matrix) for item in sublist]
    if pssm_sdt:
      pssm_sdt1 = [item for sublist in pssm_sdt_func(matrix) for item in sublist]
    
    ## add the length of matrix( protein ) as a feature
    return pssm_sdt1 + pssm_ddt1 #+ [len(matrix) ]# + trigram_1  + res_trnsfm#+ pssm_sdt1 + pssm_ddt1 #+ trigram_1 # + [random.random() for x in range(300)]



f1=open(r'GN_pssm_PA1.data','rb')
pssm1=p.load(f1)
aac=[]
RPT=[] 
#for i in range(len(pssm1)):
#    aac_pssm_obtain=aac_pssm(pssm1[i])
#    aac.append(aac_pssm_obtain)

for i in range(len(pssm1)):
    RPT_obtain=residue_probing_transform(pssm1[i])
    RPT_obtain1=np.array(RPT_obtain)
    RPT_obtain2=RPT_obtain1.T
    RPT_obtain=np.reshape(RPT_obtain2,(1,400)) 
#    for j in range(20):
#        sdt_pssm_obtain=np.insert(sdt_pssm_obtain2[j+1], 0, values=sdt_pssm_obtain2[j])
    RPT.append(RPT_obtain)
RPT1=np.array(RPT)
#sio.savemat('186_pssm_sdt.mat',{'PDB186_pssm_sdt':sdt1})
[m,n,q]=np.shape(RPT1)
X1=np.reshape(RPT1,(m,400))
sio.savemat('GN_RPT.mat',{'PDBGN_RPT':X1})












