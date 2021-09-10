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

EDT=[] 

f1=open(r'GN_pssm_PA1.data','rb')
pssm1=p.load(f1)
aac=[]
EDT=[]
#for i in range(len(pssm1)):
#    aac_pssm_obtain=aac_pssm(pssm1[i])
#    aac.append(aac_pssm_obtain)

for i in range(len(pssm1)):
    EDT_obtain=evolutionary_diff(pssm1[i])
    EDT_obtain1=np.array(EDT_obtain)
    EDT_obtain2=EDT_obtain1.T
    EDT_obtain=np.reshape(EDT_obtain2,(1,400))
    
#    for j in range(20):
#        sdt_pssm_obtain=np.insert(sdt_pssm_obtain2[j+1], 0, values=sdt_pssm_obtain2[j])
        
        

    EDT.append(EDT_obtain)

EDT1=np.array(EDT)
#sio.savemat('186_pssm_sdt.mat',{'PDB186_pssm_sdt':sdt1})
[m,n,q]=np.shape(EDT1)
X1=np.reshape(EDT1,(m,400))
sio.savemat('GN_EDT.mat',{'PDBGN_EDT':X1})












#import random
#from helpers import *
#
#import math
#import os
#
#def load_data(positive, negative):
#    random.seed(9999) #in case load data called many times
#    #random.seed(1000000002)
#    #random.seed(100)
#    #random.seed(10000)
#    #dataset to check against
#    file_names = []
#    #with open("id_list_blast.txt") as f:
#    #    file_names = [ x.split(" ")[0] for x in f ]
#    #    #print(file_names)
#    dataset = []
#    positives = []
#    negatives = []
#    #not_having_disorder = []
#    with open(positive) as file1:
#        content2 = file1.read().replace("\n", " ").split(">")
#        content = [a.split(" ", 1) for a in content2]
#        #print(len(content ) )
#        for a in content[1:]:
#          if a != "" and a!= [""] and a != []:
#            if not  os.path.exists('./downloaded_dataset/mgbg/'
#            + a[0].split("|")[0] + ".monogram" ):
#              print(a)
#            
#            if (os.path.exists('./downloaded_dataset/PSSM/'+
#            a[0].split("|")[0] + ".mat") and
#            #os.path.exists('./downloaded_dataset/DisPredict_V2.0/Software/Output/prediction/' + a[0].split("|")[0] + "/" + a[0].split("|")[0] + ".dispredict2.sl477.predict") and
#            os.path.exists('./downloaded_dataset/mgbg/' + a[0].split("|")[0] + ".monogram" ) and
#            os.path.exists('./downloaded_dataset/mgbg/' + a[0].split("|")[0] + ".bigram") ):
#                #if a[0].split("|")[0] in file_names:
#              if len(a[1].replace(" ", "") ) > 8:
#                item = a[1].replace(" ", "")
#                if 'X' not in item and 'U' not in item :#and 'O' not in item:
#                  positives.append([a[0], a[1].replace(" ", ""), 1])
#                #if not os.path.exists('./downloaded_dataset/DisPredict_V2.0/Software/Output/prediction/' + a[0].split("|")[0] + "/" + a[0].split("|")[0] + ".dispredict2.sl477.predict") :
#                #    not_having_disorder.append([a[0], a[1].replace(" ", ""),])
#    #print(not_having_disorder)
#    #print(len(not_having_disorder))
#    #with open("downloaded_dataset/remaining_binding_17_final.txt", "w") as f:
#    #    for x in not_having_disorder:
#    #        f.write(x[0] + '\n')
#    #print(len(positives) )
#    with open(negative) as file2:
#      negatives = []     
#      content2 = file2.read().replace("\n", " ").split(">")
#      
#      content = [a.split(" ", 1) for a in content2]
#      #print(len(content) ) 
#      total = len(content)
#      for a in content[1:]:
#        if a != "" and a != [""] and a != []:
#          #print(a)
#          #check if all the values have been calculated first
#          #look into all folders
#          if not os.path.exists('./downloaded_dataset/mgbg/'
#          + a[0].split("|")[0] + ".monogram" ):
#            print("hello")
#            print(a)
#          if (os.path.exists('./downloaded_dataset/PSSM/'+ a[0].split("|")[0] + ".mat") and
#          #os.path.exists('./downloaded_dataset/PSSM/'+ a[0].split("|")[0][:4]
#          #+ ":" + a[0].split("|")[0][4:]
#          #+ ".mat")  )and
#          #os.path.exists('./downloaded_dataset/DisPredict_V2.0/Software/Output/prediction/' + a[0].split("|")[0] + "/" + a[0].split("|")[0] + ".dispredict2.sl477.predict") and
#          os.path.exists('./downloaded_dataset/mgbg/' + a[0].split("|")[0] + ".monogram" ) and
#          os.path.exists('./downloaded_dataset/mgbg/' + a[0].split("|")[0]
#          + ".bigram") ) :
#            #  #and if they do not contain any illegal charecters
#            #  #print(a)
#            #  #if a[0].split("|")[0] in file_names:
#            item = a[1].replace(" ", "")
#            #print(a)
#            if len(a[1].replace(" ", "") ) > 8:
#              if 'X' not in item and 'U' not in item: # and 'O' not in item:
#                negatives.append([a[0], a[1].replace(" ", ""), 0])
#            #else:
#            #  print("hello" + str(a) )
#    print(len(negatives))
#    #print(len(negatives[:758]))
#    #now split negatives into three halves:
#    #rand_values = random.sample(range(0, len(negatives ) ), len(negatives ) )
#    #my_negatives = [negatives[x] for x in rand_values]#[800 : 2 * 800]
#    random.shuffle(negatives)
#    #random.shuffle(negatives)
#    length_negatives = 758 #len(positives) 
#    neg = negatives[: length_negatives]
#    returnVal = positives + neg #[200 : 759 ] # we have total 559 positive samples
#    random.shuffle(returnVal , random.random)
#    return returnVal
#    
#if __name__ == "__main__":
#    import helpers
#    #data = load_data("Supp-S1/binding_fasta.txt", "Supp-S1/nonbinding_fasta.txt") 
#    data = load_data("Supp-S1/binding.txt",
#    "Supp-S1/non-binding.txt") 
#    print(len([y for y in data if y[2] == 1 ]) )
#    #data = load_data("Supp-S1/all_binding_with_fasta_final.txt",
#    print(len(data) )
#    #"Supp-S1/all_non_binding_with_fasta.txt")
#    #print(len(data) )
#
#    #with open("results/data.txt", "w") as f:
#    #  for x in data:
#    #    f.write(str(x[0].split("|")[0] ) + "," + str(x[2] )  + "\n" + str(x[1].replace(" ", "") )+
#    #    "\n" )
