function [P] = execute_transform(X, Y, transform_type, transform_parameter)
%
% Input:
%      X -- a N*D training data matrix, where each row indicates an instance
%      Y -- a N*Q binary label matrix
%      transform_type -- transform type
%            =0, no transform
%            =1, PCA
%            =2, CCA
%            =3, MLSI
%            =4, MDDM
%            =5, MVMD
%            =6, MLDA
%            =7, DMLDA
%            ...
%     parameter.ratio -- a constant to detect the number of dimensions
%     parameter.beta  -- a trade-off factor in MLSI and MVMD
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameter = transform_parameter;

switch(transform_type)
    case 1 %PCA
        disp('PCA -------------------');
        [P] = PCA_transform(X,Y,parameter);
        
    case 2 %CCA
        disp('CCA--------------------');
        [P] = CCA_transform(X,Y,parameter);
 
    case 3 %MLSI
        disp('MLSI--------------------');
        [P] = MLSI_transform(X,Y,parameter);
        
    case 4 %MDDMp with orthonormal transforms
        disp('MDDMp--------------------');
        [P] = MDDM_transform(X,Y,parameter,0);
        
    case 5 %MDDMf with orthonormal features
        disp('MDDMf--------------------');
        [P] = MDDM_transform(X,Y,parameter,1);
        
    case 6 %MVMD with orthonormal transforms
        disp('MVMD--------------------');
        [P] = MVMD_transform(X,Y,parameter);   
        
    case 7 % DMLDA (Sb & Sw)
        disp('DMLDA--------------------');
        [P] = DMLDA_transform(X,Y,parameter,3);
           
    case 8  %MLDA (Sb & Sw)(wMLDAc with correlation-based weight form)
        disp('MLDA (wMLDAc) --------------------');
        [P] = wMLDA_transform(X,Y,parameter,3);
      
    case 9 %wMLDAb with binary weight form
        disp('wMLDAb--------------------');
        [P] = wMLDA_transform(X,Y,parameter,1);
        
    case 10 %wMLDAe with entropy weight form
        disp('wMLDAe--------------------');
        [P] = wMLDA_transform(X,Y,parameter,2);
        
    case 11 %wMLDAf with fuzzy weight form
        disp('wMLDAf--------------------');
        [P] = wMLDA_transform(X,Y,parameter,4);
        
    case 12 %wMLDAd with dependence weight form
        disp('wMLDAd--------------------');
        [P] = wMLDA_transform(X,Y,parameter,5);
        
end

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
