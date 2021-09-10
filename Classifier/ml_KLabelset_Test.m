 %function [exampleMetrics,labelMetrics,testLabels,labelQuantity1,labelQuantity2] = ml_KLabelset_Test( test_data,test_targets,outputStruct,KLabelsetsSelected,classLabel,threshold )

function [testLabels,labelQuantity1,labelQuantity2] = ml_KLabelset_Test( test_data,test_targets,outputStruct,KLabelsetsSelected,classLabel,threshold )

    [numSamples,numClasses]         = size(test_targets);   
    classifierNum                  	= size(outputStruct,2);
    labelQuantity1               	= zeros(numSamples,numClasses);
	labelQuantity2                	= zeros(numSamples,numClasses);
    
    % Testing
    finallabel                      = zeros(numSamples,classifierNum);
	parfor i = 1:classifierNum   
        fprintf('K-Labelset Method: testing the classifier for %d-th labelset\n',i);
      	[finallabel(:,i),~]                       	= base_svm_OAA_test( outputStruct{1,i},test_data,classLabel{1,i} );
	end
    
	for i = 1:classifierNum
       	tempLabelset                                = KLabelsetsSelected{1,i};
    	tempLabelsetBinary                          = zeros(1,numClasses);
     	tempLabelsetBinary(1,tempLabelset)          = 1;
    	labelQuantity1                              = labelQuantity1 + repmat(tempLabelsetBinary,numSamples,1);     
     	tempIndex                                   = double( dec2bin( finallabel(:,i)-1, length(tempLabelset) ) ) - '0';
      	labelQuantity2(:,tempLabelset)              = labelQuantity2(:,tempLabelset) + tempIndex;       
	end
    
	testLabels        	= ((labelQuantity2 ./ labelQuantity1) >= threshold) *2 - 1;

    % Computing the Evaluation Metrics Values 
%     [exampleMetrics,labelMetrics]   = ml_0evaluateMetrics( testLabels,test_targets,numClasses,numSamples );
    
end

