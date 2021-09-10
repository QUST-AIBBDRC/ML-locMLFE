function [outputStruct,KLabelsetsSelected,classLabel,threshold,coverage] = ml_KLabelset_active_overlap_Train( train_data,train_targets,type,ker,kerp,C,labelNum,classifierNum,beta )
% train_data: n*d, one row for an instance with d features
% train_targets: n*m, one row for an instance, -1 means irrelevant, 1 means relevant

    [~,numClasses]                  = size(train_targets);
    KLabelsets                      = nchoosek(1:numClasses,labelNum);
	classifierNum                   = min(classifierNum, size(KLabelsets,1));

    classLabel                      = cell(1,classifierNum);
    for i = 1:classifierNum
       classLabel{1,i}              = zeros(1,0); 
    end
    
    selectedIndex                   = zeros(1,classifierNum);
    numCheck                        = 5;

    outputStruct                    = cell(1,classifierNum);     
    KLabelsetsSelected              = cell(1,classifierNum); 
    labelFrequency                  = zeros(1,numClasses);
  
    addMatrix                           = (repmat((1:1:size(KLabelsets,1))',1,labelNum)-1)*numClasses;
    KLabelsetsNew                       = KLabelsets + addMatrix;
    KLabelsetsReshape                   = reshape(KLabelsetsNew',[1 size(KLabelsets,1)*labelNum]);
    KLabelsetsCode                      = zeros(1,size(KLabelsets,1)*numClasses);
    KLabelsetsCode(KLabelsetsReshape)   = 1;
    KLabelsetsCodeMatrix                = reshape(KLabelsetsCode, [numClasses size(KLabelsets,1)])';
    
    mark                                = ones(size(KLabelsets,1),1);
    
    fRatioK                             = zeros(1,size(KLabelsets,1));
    jEntropy                            = zeros(1,size(KLabelsets,1));
    quality                             = zeros(1,size(KLabelsets,1));
    coverage                            = zeros(1,classifierNum);
    
    % Training Multi-label Classifiers
    ytrain                          = zeros(size(train_data,1),classifierNum);
    for i = 1:classifierNum  
        fprintf('Overlap Active K-Labelset: training the classifier for %d-th labelset\n',i);   
        if i == 1
          	%%%%%%%%%%%%%%% Select Label Set Iteratively %%%%%%%%%%%%%
            randFirst                       = randi(size(KLabelsets,1));
            tempLabelset                    = KLabelsets(randFirst,:);
            labelFrequency(tempLabelset)    = labelFrequency(tempLabelset) + 1;
            mark(randFirst)                 = Inf;
            jEntropy(randFirst)           	= jointEntropy( train_targets, tempLabelset );
            fRatioK(randFirst)            	= fisherRatio( train_data, train_targets, tempLabelset, sqrt(1/(2*kerp)) );
            quality(randFirst)            	= beta* fRatioK(randFirst) + (1-beta) * jEntropy(randFirst);
            baseQuality                     = quality(randFirst);
            selectedIndex(1)                = randFirst;
        else
            frequencyMatrix                 = repmat(labelFrequency,size(KLabelsets,1),1);
            frequencyMatrixNew              = (frequencyMatrix + KLabelsetsCodeMatrix);
            frequencyDiff                  	= (max(frequencyMatrixNew,[],2) - min(frequencyMatrixNew,[],2)) .* mark;
            minFrequencyDiff                = min( frequencyDiff );
            bestIndices                     = find(frequencyDiff == minFrequencyDiff);
            bestNum                         = length(bestIndices);
            seedOrder                       = randperm(bestNum);
            
            for k = 1 : min(numCheck,bestNum)
         	%%%%%%%%%%%%% Calculate Joint Entropy Fisher Ratio %%%%%%%%%%%
                tempLabelset                                = KLabelsets(bestIndices(seedOrder(k)),:);
                if jEntropy(bestIndices(seedOrder(k))) == 0
                    [jEntropy(bestIndices(seedOrder(k))), jointProbsNonZeroNum]   	= jointEntropy( train_targets, tempLabelset );
                    fRatioK(bestIndices(seedOrder(k)))   	= fisherRatio( train_data, train_targets, tempLabelset, sqrt(1/(2*kerp)) );
                    quality(bestIndices(seedOrder(k)))    	= fRatioK(bestIndices(seedOrder(k))) + jEntropy(bestIndices(seedOrder(k)));
                end
                if quality(bestIndices(seedOrder(k))) > baseQuality && jointProbsNonZeroNum==2^labelNum
                 	selectedIndex(i)    = bestIndices(seedOrder(k)); 
                  	baseQuality         = max(quality(selectedIndex(1:i)));
                  	break;
                end
                if k == min(numCheck,bestNum) && selectedIndex(i) == 0
                   	[~,bestIndex]       = max( quality(bestIndices) ); 
                  	selectedIndex(i)    = bestIndices(bestIndex);
                 	tempLabelset      	= KLabelsets(selectedIndex(i),:);
                  	baseQuality         = max(quality(selectedIndex(1:i)));
                end
            end

            labelFrequency(tempLabelset)    = labelFrequency(tempLabelset) + 1;
            mark(selectedIndex(i))         	= Inf;
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
        coverage(1,i)               = length(unique(ytrain(:,i))) / 2^labelNum;
    end
    
    % Finding the best threshold
	candidates                              = [0.3 0.5 0.7]; 
	Fscores                                 = zeros(1,length(candidates));
	[~,~,~,labelQuantity1,labelQuantity2] 	= ml_KLabelset_Test( train_data,train_targets,outputStruct,KLabelsetsSelected,classLabel,0.5 );
	parfor j = 1:length(Fscores)
     	testLabels              = ((labelQuantity2 ./ labelQuantity1) >= candidates(j)) *2 - 1;
     	examPrecision           = mean(sum(testLabels == 1 & train_targets == 1, 2) ./ (sum(testLabels == 1, 2) + eps));
    	examRecall              = mean(sum(testLabels == 1 & train_targets == 1, 2) ./ (sum(train_targets == 1, 2) + eps));
       	Fscores(j)              = (2 * examPrecision * examRecall) / (examPrecision + examRecall + eps);
	end

	[~,index]                   = max(Fscores);
	threshold                   = candidates(index);

end


