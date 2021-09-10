function [outputStruct,KLabelsetsSelected,classLabel,coverage] = ml_KLabelset_active_disjoint_Train( train_data,train_targets,type,ker,kerp,C,labelNum,beta )
% train_data: n*d, one row for an instance with d features
% train_targets: n*m, one row for an instance, -1 means irrelevant, 1 means relevant

    [~,numClasses]                  = size(train_targets);
    WholeLabelsets                 	= 1:1:numClasses;
    classifierNum                   = floor(numClasses/labelNum);  

    KLabelsetsSelected              = cell(1,classifierNum);  
    classLabel                      = cell(1,classifierNum);
    for i = 1:classifierNum
       classLabel{1,i}              = zeros(1,0); 
    end

	coverage                      	= zeros(1,classifierNum);
    
    outputStruct                    = cell(1,classifierNum);   
    currentLabelsets              	= WholeLabelsets;
    % Training Multi-label Classifiers
    ytrain                          = zeros(size(train_data,1),classifierNum);
    for i = 1:classifierNum  
        fprintf('Disjoint Active K-Labelset: training the classifier for %d-th labelset\n',i);           
        if i ~= classifierNum
            randFirst               = randperm(length(currentLabelsets),1);
            tempLabelset          	= currentLabelsets(randFirst);
            currentLabelsets        = setdiff(currentLabelsets,currentLabelsets(randFirst));
            %%%%%%%%%%%%% Calculate Joint Entropy Fisher Ratio %%%%%%%%%%%
            for k = 2:labelNum 
                jEntropy            = zeros(1,length(currentLabelsets));
                fRatioK             = zeros(1,length(currentLabelsets));
                for q = 1:length(currentLabelsets)
                    labels          = [tempLabelset currentLabelsets(q)];
                    jEntropy(q)     = jointEntropy( train_targets, labels );
                    fRatioK(q)      = fisherRatio( train_data, train_targets, labels, sqrt(1/(2*kerp)) );
                end
%                 fRatioK             = mapminmax(fRatioK,0,1);
                quality             = beta * fRatioK + (1-beta) * jEntropy;
                [~,index]           = max(quality);
                tempLabelset      	= [tempLabelset currentLabelsets(index)];
                currentLabelsets   	= setdiff(currentLabelsets,currentLabelsets(index));
            end
        else
            tempLabelset           	= currentLabelsets;
        end

        train_targets_temp          = train_targets(:,tempLabelset);
        for j = 0:2^length(tempLabelset)-1
          	tempIndex               = double( dec2bin(j,length(tempLabelset)) ) - '0';
          	tempIndexMatrix         = repmat( tempIndex, size(train_data,1), 1 );
          	flags               	= sum( (train_targets_temp+1)/2 == tempIndexMatrix , 2) == length(tempLabelset);
            if sum(flags) ~= 0
              	ytrain(flags,i)    	= j+1;
             	classLabel{1,i}     = [classLabel{1,i}, j+1];
            end
        end
        KLabelsetsSelected{1,i}     = tempLabelset;
    end
    
    parfor i = 1:classifierNum
       outputStruct{1,i}           = base_svm_OAA_train( [train_data ytrain(:,i)],type,ker,kerp,C ); 
       coverage(1,i)               = length(unique(ytrain(:,i))) / 2^(length(KLabelsetsSelected{1,i}));
    end
    
end


