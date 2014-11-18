function [countsMatrix wrongPredictionsNodes wrongPredictionsWeights wrongPredictionsCorrect] = computeStats(matrix1,matrix2,matrix3)
    possibilities=cell(3,3);
    possibilities{1,1}='True Positive';
    possibilities{1,2}='False Negative';
    possibilities{1,3}='False Positive';
    possibilities{2,1}='False Positive';
    possibilities{2,2}='True Negative';
    possibilities{2,3}='False Positive';
    possibilities{3,1}='False Positive';
    possibilities{3,2}='False Negative';
    possibilities{3,3}='True Positive';
    countsMatrix=zeros(3,3);
    wrongPredictionsNodes={};
    wrongPredictionsWeights=[];
    wrongPredictionsCorrect=[];
    for i=1:length(nodeNames)
        for j=1:length(nodeNames)
            ijthPossibility=possibilities{matrix1(i,j)+2,matrix2(i,j)+2};
            countsMatrix(matrix1(i,j)+2,matrix2(i,j)+2)= ...
                countsMatrix(matrix1(i,j)+2,matrix2(i,j)+2)+1;
            if(strcmp(ijthPossibility,'False Positive') || strcmp(ijthPossibility,'False Negative'))
                ijthPossibility=possibilities{matrix3(i,j)+2,matrix2(i,j)+2};
                countsMatrix(matrix1(i,j)+2,matrix2(i,j)+2)= ...
                    countsMatrix(matrix1(i,j)+2,matrix2(i,j)+2)-1;
                countsMatrix(matrix3(i,j)+2,matrix2(i,j)+2)= ...
                    countsMatrix(matrix3(i,j)+2,matrix2(i,j)+2)+1;
                if(strcmp(ijthPossibility,'False Positive') || strcmp(ijthPossibility,'False Negative'))
                    wrongPredictionsNodes{end+1}=[metsNames{i} ' ' metsNames{j}];
                    wrongPredictionsWeights(end+1)=matrix2Continuous(i,j);
                    wrongPredictionsCorrect(end+1)=matrix3(i,j);
                end
            end
        end
    end
end