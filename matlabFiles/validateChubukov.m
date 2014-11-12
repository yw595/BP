function []=validateChubukov(filterUnclear,dataFile,isMets)

    [junk1 junk2 totalData]=xlsread('../Chubukov/inline-supplementary-material-2.xlsx',4);

    [metsNames rxnNames]=renameMetsAndRxns();
    if(filterUnclear)
        rxnNames=rxnNames(9:end-12);
    end
    
    %parse equations written out in NodeEqs.txt
    SMatrix=formSMatrix(metsNames);

    %form NodeInfluenceMatrix by dot-product of S-matrix columns for reaction,
    %filter out elements at CO2 rows, and remove Fba and Tpi interactions since
    %are ambiguous
    rxnInfluenceMatrix=formInfluenceMatrix(fluxNames,SMatrix,metsNames);
    
    metsAndRxnsInfluenceMatrix=zeros(length(metsNames)+length(fluxNames),length(metsNames)+length(fluxNames));
    metsAndRxnsInfluenceMatrix(1:length(metsNames),length(metsNames)+1:length(metsNames)+length(fluxNames))=...
    -sign(SMatrix);
    metsAndRxnsInfluenceMatrix=metsAndRxnsInfluenceMatrix+metsAndRxnsInfluenceMatrix';
    
    nodeInfluenceMatrix=[];
    nodeNames={};
    if(isMets)
        nodeInfluenceMatrix=rxnInfluenceMatrix;
        nodeNames=[metsNames rxnNames];
    else
        nodeInfluenceMatrix=metsAndRxnsInfluenceMatrix;
        nodeNames=rxnNames;
    end

    storeNodeInfluenceMatrix=nodeInfluenceMatrix;
    absNodeInfluenceMatrix=abs(nodeInfluenceMatrix);

    %from edges file, read in sign of influence in InfluenceMatrix, size of
    %influence in InfluenceMatrixContinuous
    WNodeInfluenceMatrix=zeros(size(nodeInfluenceMatrix,1),size(nodeInfluenceMatrix,2));
    WNodeInfluenceMatrixContinuous=zeros(size(nodeInfluenceMatrix,1),size(nodeInfluenceMatrix,2));
    WFID=fopen(dataFile);
    line=fgetl(WFID);
    while(line~=-1)
        firstNode=line(1:regexp(line,'(')-2);
        secondNode=line(regexp(line,')')+2:regexp(line,'=')-2);
        influence=line(regexp(line,'(')+1:regexp(line,')')-1);
        if(sum(strcmp(secondNode,nodeNames))~=0 && sum(strcmp(firstNode,nodeNames))~=0)
            WNodeInfluenceMatrixContinuous(strcmp(nodeNames,firstNode),strcmp(nodeNames,secondNode))=...
            str2num(line(min(regexp(line,'[-|\d]+\.\d+')):end));
            if(strcmp(influence,'activates'))
                WNodeInfluenceMatrix(strcmp(nodeNames,firstNode),strcmp(nodeNames,secondNode))=1;
            elseif(strcmp(influence,'inhibits'))
                WNodeInfluenceMatrix(strcmp(nodeNames,firstNode),strcmp(nodeNames,secondNode))=-1;
            end
        end
        line=fgetl(WFID);
    end

    for z=1:10 
        %form NodeInfluenceMatrix with all influences up to distance of cutoff
        %allowed
        cutoff=z;
        for i=1:length(nodeNames)
            for j=1:length(nodeNames)
                [junk1,path,junk2]=graphshortestpath(sparse(absNodeInfluenceMatrix),i,j);
                if(length(path)<=cutoff+1 && length(path)>1)
                    sign=1;
                    for k=1:length(path)-1
                        sign=sign*storeNodeInfluenceMatrix(path(k),path(k+1));
                    end
                    nodeInfluenceMatrix(i,j)=sign;
                end
            end
        end
    
        %iterate over NodeInfluenceMatrix and W, calculate accuracy statistics
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
                ijthPossibility=possibilities{storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2};
                countsMatrix(storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)= ...
                countsMatrix(storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)+1;
                %if false positive or negative, use extended NodeInfluence instead,
                %adjusting countsMatrix. If still false, write out Nodes, weight of
                %predicted interaction, and correct interaction
                if(strcmp(ijthPossibility,'False Positive') || strcmp(ijthPossibility,'False Negative'))
                    ijthPossibility=possibilities{NodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2};
                    countsMatrix(storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)= ...
                    countsMatrix(storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)-1;
                    countsMatrix(nodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)= ...
                    countsMatrix(nodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)+1;
                    if(strcmp(ijthPossibility,'False Positive') || strcmp(ijthPossibility,'False Negative'))
                        wrongPredictionsNodes{end+1}=[metsNames{i} ' ' metsNames{j}];
                        wrongPredictionsWeights(end+1)=WNodeInfluenceMatrixContinuous(i,j);
                        wrongPredictionsCorrect(end+1)=nodeInfluenceMatrix(i,j);
                    end
                end
            end
        end
    
        %display wrong Predictions
        for i=1:length(wrongPredictionsNodes)
            if(wrongPredictionsWeights(i)~=0 && wrongPredictionsCorrect(i)~=0)
                disp([wrongPredictions{i} ' ' num2str(wrongPredictionsWeights(i)) ...
                ' ' num2str(wrongPredictionsCorrect(i))]);
            end
        end
    
        recall=(countsMatrix(1,1)+countsMatrix(3,3))/(sum(countsMatrix(1,:)+countsMatrix(3,:)));
        precision=(countsMatrix(1,1)+countsMatrix(3,3))/(sum(countsMatrix(:,1)+countsMatrix(:,3)));
        disp(z)
        disp(recall)
        disp(precision)
        countsMatrix
    end
end