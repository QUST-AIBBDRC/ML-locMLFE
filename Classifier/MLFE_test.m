function  [Outputs,Pre_Labels] = MLFE_test(testX, trainX, test_target, model)

%MLFE_TEST 
%    Syntax
%
%      [OneError,Coverage,RankingLoss,Average_Precision,MacF1,MicF1,Outputs,Pre_Labels] = MLFE_test(testX, trainX,test_target, model)
%
  %  Description
%
%           MLFE_TEST  takes,
         %  testX --- An P1xD array, the ith instance of testing instance is stored in train_data(i,:)
%           trainX--- An P2xD array, the ith instance of training instance is stored in train_data(i,:)
          % test_target --- A P1xQ array, if the ith training instance belongs to the jth class, then train_target(i,j) equals +1, otherwise train_target(i,j) equals -1
           %model   --- trained model parameters
% 
%       and returns,
%           RankingLoss      - The ranking loss on testing data
%           OneError         - The one-error on testing data
%           Coverage         - The coverage on testing data
%           Average_Precision- The average precision on testing data
%           MacF1            - The macro-averaging f1 on testing data
%           MicF1            - The micro-averaging f1 on testing data
%           Outputs          - A QxP1 array, the probability of the ith testing instance belonging to the jCth class is stored in Outputs(j,i)
%           Pre_Labels       - A QxP1 array, if the ith testing instance belongs to the jth class, then Pre_Labels(j,i) is +1, otherwise Pre_Labels(j,i) is -1

% Copyright: Qian-Wen Zhang (cowenzhang@tencent.com,zhangqw@seu.edu.cn),
%   Yun Zhong(zeuszhong@tencent.com),
%   Min-Ling Zhang (mlzhang@seu.edu.cn)
%   Tencent SPPD, Chengdu 610041, China  &&
%   School of Computer Science and Engineering, Southeast University,
%   Nanjing 211189, China
%



fprintf(1,'\nPredict multi-label for the test data.\n');
%Compute kernel matrix for prediction using testX and trainX
Ktest = kernelmatrix(model.ker, model.par, testX, trainX);
%Prediction.
degree = Ktest*model.Beta+repmat(model.b,size(Ktest,1),1);

label = zeros(size(degree));
label(degree >= 0) = 1;
label(degree < 0) = -1;
Pre_Labels=label';

Max=max(max(degree));
Min=min(min(degree));
Outputs=(degree-Min)/(Max-Min);
Outputs=Outputs';
% RankingLoss=Ranking_loss(Outputs,test_target);
% OneError=One_error(Outputs,test_target);
% Coverage=coverage(Outputs,test_target);
% Average_Precision=Average_precision(Outputs,test_target);
% MacF1 = MacroF1(Pre_Labels,test_target);
% MicF1 = MicroF1(Pre_Labels,test_target);

end