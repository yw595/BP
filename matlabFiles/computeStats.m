function [countsMatrix wrongPredictionsNodes wrongPredictionsWeights wrongPredictionsCorrect] = ...
    computeStats(matrix1,matrix2,matrix3,unsigned)

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
    if(exist('unsigned','var') && unsigned)
        possibilities=cell(2,2);
        possibilities{1,1}='True Positive';
        possibilities{1,2}='False Negative';
        possibilities{2,1}='False Positive';
        possibilities{2,2}='True Negative';
    end
    
    countsMatrix=zeros(3,3);
    wrongPredictionsNodes={};
    wrongPredictionsWeights=[];
    wrongPredictionsCorrect=[];
    for i=1:length(nodeNames)
        for j=1:length(nodeNames)
            
            coordinate1 = matrix1(i,j)+2;
            coordinate2 = matrix2(i,j)+2;
            if(exist('unsigned','var') && unsigned)
                coordinate1 = abs(matrix1(i,j));
                if(coordinate1==0)
                    coordinate1 = 2;
                end
                coordinate2 = abs(matrix2(i,j));
                if(coordinate2==0)
                    coordinate2 = 2;
                end
            end
            
            ijthPossibility=possibilities{coordinate1,coordinate2};
            countsMatrix(coordinate1,coordinate2)= ...
                countsMatrix(coordinate1,coordinate2)+1;
            if(strcmp(ijthPossibility,'False Positive') || strcmp(ijthPossibility,'False Negative'))
                ijthPossibility=possibilities{matrix3(i,j)+2,coordinate2};
                countsMatrix(coordinate1,coordinate2)= ...
                    countsMatrix(coordinate1,coordinate2)-1;
                
                coordinate3 = matrix3(i,j)+2;
                coordinate2 = matrix2(i,j)+2;
                if(exist('unsigned','var') && unsigned)
                    coordinate3=abs(matrix1(i,j));
                    if(coordinate3==0)
                        coordinate3 = 2;
                    end
                    coordinate2=abs(matrix2(i,j));
                    if(coordinate2==0)
                        coordinate2 = 2;
                    end
                end
                
                countsMatrix(coordinate3,coordinate2)= ...
                    countsMatrix(coordinate3,coordinate2)+1;
                if(strcmp(ijthPossibility,'False Positive') || strcmp(ijthPossibility,'False Negative'))
                    wrongPredictionsNodes{end+1}=[metsNames{i} ' ' metsNames{j}];
                    wrongPredictionsWeights(end+1)=matrix2Continuous(i,j);
                    wrongPredictionsCorrect(end+1)=matrix3(i,j);
                end
            end
        end
    end
end