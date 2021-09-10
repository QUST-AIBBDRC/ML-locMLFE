function [Outputs,Pre_Labels]=InsDif(train_data,train_target,test_data,test_target,ratio)
%InsDif tackles multi-label learning problem through instance differentiation
%
%    Syntax
%
%       [HammingLoss,RankingLoss,OneError,Coverage,Average_Precision,Outputs,Pre_Labels,num_iter]=InsDif(train_data,train_target,test_data,test_target,ratio)
%
%    Description
%
%       MIML takes,
%           train_data       - An M1xN array, the ith instance of training instance is stored in train_data(i,:)
%           train_target     - A QxM1 array, if the ith training instance belongs to the jth class, then train_target(j,i) equals +1, otherwise train_target(j,i) equals -1
%           test_data        - An M2xN array, the ith instance of testing instance is stored in test_data(i,:)
%           test_target      - A QxM2 array, if the ith testing instance belongs to the jth class, test_target(j,i) equals +1, otherwise test_target(j,i) equals -1
%           ratio            - The number of clusters is set to ratio*M1

%      and returns,
%           HammingLoss      - The hamming loss on testing data as described in [1]
%           RankingLoss      - The ranking loss on testing data as described in [1]
%           OneError         - The one-error on testing data as described in [1]
%           Coverage         - The coverage on testing data as described in [1]
%           Average_Precision- The average precision on testing data as described in [1]
%           Outputs          - The output of the ith testing instance on the jth class is stored in Outputs(j,i)
%           Pre_Labels       - If the ith testing instance belongs to the jth class, then Pre_Labels(j,i) is +1, otherwise Pre_Labels(j,i) is -1
%           num_iter         - Number of iterations in first layer clustering
%
% [1] M.-L. Zhang and Z.-H. Zhou. Multi-label learning by instance differentiation. In: Proceedings of the 22nd AAAI Conference on Artificial Intelligence (AAAI'07), Vancouver, Canada, 2007.
% [2] Schapire R. E., Singer Y. BoosTexter: a boosting based system for text categorization. Machine Learning, 39(2-3): 135-168, 2000.
%
     
     %Preparing data
     [num_class,num_train]=size(train_target);
     [num_test,Dim]=size(test_data);
     Outputs=zeros(num_class,num_test);
     Pre_Labels=zeros(num_class,num_test);
     num_cluster=floor(ratio*num_train);
     
     t_vec=rand(num_class,Dim);
     for iter=1:num_class
        index=[];
        for i=1:num_train
            if(train_target(iter,i)==1)
                index=[index,i];
            end
        end
        if(~isempty(index))
            all=train_data(index,:);
            t_vec(iter,:)=mean(all,1);
        end
    end
     
     %training     
     
     inner_train_train=zeros(num_train,num_train);
     inner_test_train=zeros(num_test,num_train);
     inner_train_center=zeros(num_train,num_class);
     inner_test_center=zeros(num_test,num_class);
     inner_center_center=zeros(num_class,num_class);
     norm_train=zeros(1,num_train);
     norm_test=zeros(1,num_test);
     num_center=zeros(1,num_class);

     norm_train=sum((train_data.*train_data)');
     norm_test=sum((test_data.*test_data)');
     for i=1:num_class
         norm_center(1,i)=t_vec(i,:)*t_vec(i,:)';
     end

     disp('Computing matrix inner_train_train');
     for i=1:num_train
         if(mod(i,100)==0)
             disp(strcat('Computing inner_train_train for train instance:',num2str(i)));
         end
         inner_train_train(i,i)=norm_train(1,i);
         for j=(i+1):num_train
             inner_train_train(i,j)=train_data(i,:)*train_data(j,:)';
             inner_train_train(j,i)=inner_train_train(i,j);
         end
     end

     disp('Computing matrix inner_test_train');
     for i=1:num_test
         if(mod(i,100)==0)
             disp(strcat('Computing inner_test_train for test instance:',num2str(i)));
         end
         for j=1:num_train
             inner_test_train(i,j)=test_data(i,:)*train_data(j,:)';
         end
     end

     disp('Computing matrix inner_train_center');
     for i=1:num_train
         for j=1:num_class
             inner_train_center(i,j)=train_data(i,:)*t_vec(j,:)';
         end
     end

     disp('Computing matrix inner_test_center');
     for i=1:num_test
         for j=1:num_class
             inner_test_center(i,j)=test_data(i,:)*t_vec(j,:)';
         end
     end

     disp('Computing matrix inner_center_center');
     for i=1:num_class
         for j=1:num_class
             inner_center_center(i,j)=t_vec(i,:)*t_vec(j,:)';
         end
     end

     disp('Computing distance_matrix for clustering');


     distance_matrix=zeros(num_train,num_train);
     for i=1:(num_train-1)
         if(mod(i,100)==0)
             disp(strcat('Computing distance for train bags:',num2str(i)));
         end
         for j=(i+1):num_train
             dist=zeros(num_class,num_class);
             for m=1:num_class
                 for n=1:num_class
                     temp1=norm_train(1,i)+norm_center(1,m)-2*inner_train_center(i,m);
                     temp2=norm_train(1,j)+norm_center(1,n)-2*inner_train_center(j,n);
                     temp3=-2*(inner_train_train(i,j)-inner_train_center(i,n)-inner_train_center(j,m)+inner_center_center(m,n));
                     dist(m,n)=sqrt(temp1+temp2+temp3);
                 end
             end
             distance_matrix(i,j)=max(max(min(dist)),max(min(dist')));
         end
     end
     distance_matrix=distance_matrix+distance_matrix';

     [clustering,matrix_fai,num_iter]=MIML_cluster(num_cluster,distance_matrix);      
     Weights=real(pinv(matrix_fai)*(train_target'));
     
     %testing
     for i=1:num_test         
         if(mod(i,100)==0)
             disp(strcat('Testing for instance:',num2str(i)));
         end
         tempvec=zeros(1,num_cluster);
         for j=1:num_cluster
             index=clustering{j,1};
             dist=zeros(num_class,num_class);
             for m=1:num_class
                 for n=1:num_class
                     temp1=norm_test(1,i)+norm_center(1,m)-2*inner_test_center(i,m);
                     temp2=norm_train(1,index)+norm_center(1,n)-2*inner_train_center(index,n);
                     temp3=-2*(inner_test_train(i,index)-inner_test_center(i,n)-inner_train_center(index,m)+inner_center_center(m,n));
                     dist(m,n)=sqrt(temp1+temp2+temp3);
                 end
             end
             tempvec(1,j)=max(max(min(dist)),max(min(dist')));
         end
         Outputs(:,i)=(tempvec*Weights)';
     end
     
     for i=1:num_test
         for j=1:num_class
             if(Outputs(j,i)>=0)
                 Pre_Labels(j,i)=1;
             else
                 Pre_Labels(j,i)=-1;
             end
         end
     end

     HammingLoss=Hamming_loss(Pre_Labels,test_target);
    
     RankingLoss=Ranking_loss(Outputs,test_target);
     OneError=One_error(Outputs,test_target);
     Coverage=coverage(Outputs,test_target);
     Average_Precision=Average_precision(Outputs,test_target);