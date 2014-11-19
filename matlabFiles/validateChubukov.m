function []=validateChubukov(filterUnclear,dataFile,isMets)

    [junk1 junk2 totalData]=xlsread('../Chubukov/inline-supplementary-material-2.xlsx',4);

    [metsNames rxnNames]=renameMetsAndRxns();
    if(filterUnclear)
        rxnNames=rxnNames(9:end-12);
    end
    
    SMatrix=formSMatrix(metsNames);

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

    WNodeInfluenceMatrixContinuous=readMatrix('W_average.txt','  ');
    WNodeInfluenceMatrix=sign(WNodeInfluenceMatrixContinuous);
    
    extendedInfluenceMatrices = makeExtendedInfluenceMatrices(originalInfluenceMatrix,10);

    for z=1:10 
        
        nodeInfluenceMatrix = extendedInfluenceMatrices{z};
    
        [countsMatrix wrongPredictionsNodes wrongPredictionsWeights wrongPredictionsCorrect] = ...
            computeStats(WNodeInfluenceMatrix,storeNodeInfluenceMatrix,nodeInfluenceMatrix);
    
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