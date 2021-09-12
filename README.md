# Intro  
ML-locMLFE is a novel prediction multi-label protein subcellular localization  model, which is able to use sparse reconstruction of the training samples to represent the bottom layer of the feature space. [ML-locMLFE_framework](https://github.com/QUST-AIBBDRC/ML-locMLFE/blob/master/ML-locMLFE_framework.png)  

# System requirement  
ML-locMLFE is developed under Windows environment with:  
matlab 2015b
python  3.6.1  
numpy  1.16.1  
pandas  0.20.1   

# Build dataset
1.Build newPlant dataset following (https://www.uniprot.org/uniprot/)
2.Build newPlant dataset following (https://www.covid19dataportal.org/)

# Dataset, feature and model
The dataset file contains Gram-positive bacteria dataset, newplant dataset, SARS-CoV-2 dataset, virus dataset and Gram-negative bacteria dataset.
Feature extraction: mainEBGW.m, ebgw1.m, ebgw2.m, ebgw3.m  is the implementation of EBGW. PAAC.m,mainpseaac.m is the implementation of PseAAC. EDT.py is the implementation of EDT. RPT.py is the implementation of RPT. MCD.m, MCD1D.m, MCD2D.m, MCD3D.m, MCDexchange.m, MCDfeature.m, MCDfen.m, MCDtransform.m, MCDZD.m is the implementation of MCD. Gene Ontology can be found from http://www.ebi.ac.uk/GOA/.
Dimensional reduction: MLSI_transform.m represents MLSI. MDDM_transform.m represents MDDM. MVMD_transform.m represents MVMD. PCA_transform.m represents PCA. GRRO.m represents GRRO. MDFS.m represents MDFS.
Classifier: MLFE_train.m, MLFE_test is the implementation of MLFE. LIFT.m is the implementation of LIFT. MLKNN_train.m, MLKNN_test.m are the implementation of MLKNN. ML_RBF_train.m, ML_RBF_test.m is the implementation of ML_RBF. RankSVM_train.m, RankSVM_test.m is the implementation of RankSVM. InsDif.m is the implementation of InsDif.
Model: An example is included in the Demo file, and you can run the demo.m in MATLAB.

# contact   
Contact: 
Bin Yu  (yubin@qust.edu.cn)


