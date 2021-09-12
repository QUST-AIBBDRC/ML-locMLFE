clc,clear all

transform_type=3; %降维方法
transform_parameter.ratio =80; %降的维数
%% 融合
load('goxinxi.mat');
   shuju=[gozhengli];
%load('38testlabel.mat');
load('goxinxi.mat');

L=nellabel;
[YANGBEN,WEISHU]=size(L);
%X=[shuju,jieguo,gozhengli];%3种算法的融合
% shu=zscore(shushu);
%数据处理
train_data = shuju;
train_label = L;
[N, D] = size(train_data);
mean_train = mean(train_data);
train_CX = train_data - repmat(mean_train,N,1); %centralizing training features
mean_label = mean(train_label);
train_CY = train_label - repmat(mean_label, N, 1);
currentTrainData = train_CX;
%降维方法
transform_parameter.rank = 1;
switch(transform_type)
    case 1 %PCA

    case 2 %CCA
        %KBS
        transform_parameter.regY = 0.1;
        transform_parameter.regX = 0.1;
        
        %Neurocomputing
        %transform_parameter.regY = 0.0;
        %transform_parameter.regX = 0.0;
        
    case 3 %MLSI
        transform_parameter.beta = 0.5;
        transform_parameter.regXY = 0.1;
        transform_parameter.regX  = 0.1;
        
    case 4 %MDDMp
        
    case 5 %MDDMf
        %KBS
        transform_parameter.beta = 0.5;
        %Neurocomputing 
        %transform_parameter.beta = 1.0;
        
    case 6 %MVMD
        transform_parameter.beta  = 0.5;
        
    otherwise % wMLDA
        %KBS
        transform_parameter.regXY = 0.1;
        
        %Neurocomputing
        %transform_parameter.regXY = 0.0;
%  transform_cell{1,7}  = '-DMLDA';
% transform_cell{1,8}  = '-wMLDAc(MLDA)'; %Wang2010
% 
% transform_cell{1,9}  = '-wMLDAb';       %Park2007
% transform_cell{1,10} = '-wMLDAe';       %Xu2017 & Chen2007
% transform_cell{1,11} = '-wMLDAf';       %Xu2017 & Lin2013
% transform_cell{1,12} = '-wMLDAd';       %Xu2017
end
%%%%%%%%降维默认参数
% transform_parameter.rank =0; % 
% transform_parameter.regX = 0.1; %
% transform_parameter.regY = 0.1; %
% transform_parameter.regXY = 0.1; %
% transform_parameter.beta = 0.5; % 
%降维操作
if (transform_type >0)
    % execute a FE method
    if(transform_type >=1 & transform_type <=6) %PCA, CCA, MLSI, MDDMp/f, MVMD which need centered training labels
            [PPP] = execute_transform(currentTrainData, train_CY, transform_type, transform_parameter); 
    else % LDA type transforms which do not need centered training labels
            [PPP] = execute_transform(currentTrainData, train_label, transform_type, transform_parameter);
    end
    % convert training and testing features
    current = currentTrainData * PPP;        
end

%% 杰克刀测试+ML-kNN
%Set the ratio parameter used by LIFT
ratio=0.1;Q=[];
% Set the kernel type used by Libsvm
% svm.type='RBF';
% svm.para=[0.1];
P=[];
para.tol  = 1e-5; %tolerance
para.epsi = 0.001; %threshold paramenter, more than epsi should be penalized
para.beta1    =10;  %penalty parameter                  
para.beta2    =15; 
para.beta3    =10;
para.ker  = 'rbf'; %type of kernel function


C.c1=1; %control parameters
C.c2=2;
%  Num=9;
% Smooth=1;%KNN参数

%Smooth=1;
for i=1:519
    A=current;
    B=L';
    test_data=A(i,:);test_target=B(:,i);
    A(i,:)=[];B(:,i)=[];
    train_data=A;train_target=B;
    para.par  =1*mean(pdist(train_data));
 % [Prior,PriorN,Cond,CondN]=MLKNN_train(train_data,train_target,Num,Smooth);
  %[HammingLoss,RankingLoss,OneError,Coverage,Average_Precision,Outputs,Pre_Labels]=MLKNN_test(train_data,train_target,test_data,test_target,Num,Prior,PriorN,Cond,CondN);
%[Outputs,Pre_Labels]=MLKNN_test(train_data,train_target,test_data,test_target,Num,Prior,PriorN,Cond,CondN);
  %     clear A B train_data train_target test_data test_targe
%    [nets,errors]=BPMLL_train(train_data,train_target,100);
%    [HammingLoss,RankingLoss,OneError,Coverage,Average_Precision,Outputs,Threshold,Pre_Labels]=BPMLL_test(train_data,train_target,test_data,test_target,nets{100,1});
%[HammingLoss,RankingLoss,OneError,Coverage,Average_Precision,Outputs,Pre_Labels]=LIFT(train_data,train_target,test_data,test_target,ratio,svm);
%P=[P;Pre_Labels'];
 model = MLFE_train(train_data,train_target',para,C);
    %[OneError,Coverage,RankingLoss,Average_Precision,MacF1,MicF1,Outputs,Pre_Labels] = MLFE_test(testX, trainX,test_target, model);
[Outputs,Pre_Labels] = MLFE_test(test_data, train_data, test_target, model);
P=[P;Pre_Labels'];
Q=[Q;Outputs'];
clear A B
    
end

%% 评价指标
HL=Hamming_loss(P',L')
AP=Average_precision(Q',L')
CV=coverage(Q',L')
RL=Ranking_loss(Q',L')
OAA=0;
for i=1:519
if P(i,:)==L(i,:);
OAA=OAA+1;
end
end
zuiOAA=OAA;
%clear OAA
zuiOAA=zuiOAA/519;

