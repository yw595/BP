function [signedRecalls signedPrecisions signedCountsMatrices unsignedRecalls unsignedPrecisions unsignedCountsMatrices]= ...
    validateCorrelate(dataFile,idxsToSkip1,idxsToSkip2,isMets)
    corrThresh=0.5;

    %get reaction fluxNames in storefluxNames, take out spaces for fluxNames
    [metsNames fluxNames]=renameMetsAndRxns();
    fluxNamesIdxs=1:length(fluxNames);
    fluxNamesIdxs(idxsToSkip2)=[];
    fluxNames=fluxNames(fluxNamesIdxs);

    [fluxData fluxStds metsData metsStds] = readFluxAndMetsData(idxsToSkip1,idxsToSkip2);

    fluxCorrMatrix = makeCorrMatrix(fluxData);
    fluxCorrMatrix( fluxCorrMatrix<=0.7 ) = 0;
    
    fluxMetsCorrMatrix = makeCorrMatrix(metsData,fluxData);
    
    [junk sortIdxs] = sort(abs( reshape(fluxMetsCorrMatrix,[numel(fluxMetsCorrMatrix) 1]) ));
    fluxMetsCorrsMatrix = fluxMetsCorrsMatrix(sortIdxs(end-99:end),:);
    
    if(isMets)
        tempNodeInfluenceMatrix=sign(fluxMetsCorrMatrix);
        nodeInfluenceMatrix=zeros(length(metsNames)+length(fluxNames),length(metsNames)+length(fluxNames));
        nodeInfluenceMatrix(1:length(metsNames),length(metsNames)+1:length(metsNames)+length(fluxNames))=...
        tempNodeInfluenceMatrix;
        nodeInfluenceMatrix=nodeInfluenceMatrix+nodeInfluenceMatrix';
        nodeNames=[metsNames fluxNames];
    else
        nodeInfluenceMatrix=sign(fluxCorrMatrix);
        nodeNames=fluxNames;
    end
    
    storeNodeInfluenceMatrix=nodeInfluenceMatrix;
    absNodeInfluenceMatrix=abs(nodeInfluenceMatrix);

    WNodeInfluenceMatrixContinuous=readMatrix('W_average.txt','  ');
    WNodeInfluenceMatrix=sign(WNodeInfluenceMatrixContinuous);
    
    if(isMets)
        writeMatrix(WNodeInfluenceMatrix (WNodeInfluenceMatrix ~=nodeInfluenceMatrix),'diffMatrix.sif',0,1,0,0,metsNames,fluxNames);
        writeMatrix(WNodeInfluenceMatrix (WNodeInfluenceMatrix ~=nodeInfluenceMatrix),'diffMatrixEdges.txt',0,0,1,0,metsNames,fluxNames);
    else
        writeMatrix(WNodeInfluenceMatrix (WNodeInfluenceMatrix ~=nodeInfluenceMatrix),'diffMatrix.sif',0,1,0,0,fluxNames,fluxNames);
        writeMatrix(WNodeInfluenceMatrix (WNodeInfluenceMatrix ~=nodeInfluenceMatrix),'diffMatrixEdges.txt',0,0,1,0,fluxNames,fluxNames);
    end
    
    extendedInfluenceMatrices = makeExtendedInfluenceMatrices(originalInfluenceMatrix,10);
    
    signedRecalls=[];signedPrecisions=[];unsignedRecalls=[];unsignedPrecisions=[];signedCountsMatrices={};unsignedCountsMatrices={};
    for z=1:10 
        nodeInfluenceMatrix = extendedInfluenceMatrices{z};
    
        [countsMatrix wrongPredictionsNodes wrongPredictionsWeights wrongPredictionsCorrect] = ...
            computeStats(WNodeInfluenceMatrix,storeNodeInfluenceMatrix,nodeInfluenceMatrix);
    
        recall=(countsMatrix(1,1)+countsMatrix(3,3))/(sum(countsMatrix(1,:)+countsMatrix(3,:)));
        precision=(countsMatrix(1,1)+countsMatrix(3,3))/(sum(countsMatrix(:,1)+countsMatrix(:,3)));
        signedRecalls(z)=recall;
        signedPrecisions(z)=precision;
        signedCountsMatrices{z}=countsMatrix;
        
        [countsMatrix wrongPredictionsNodes wrongPredictionsWeights wrongPredictionsCorrect] = ...
            computeStats(WNodeInfluenceMatrix,storeNodeInfluenceMatrix,nodeInfluenceMatrix,1);
    
        recall=countsMatrix(1,1)/sum(countsMatrix(1,:));
        precision=countsMatrix(1,1)/sum(countsMatrix(:,1));
        unsignedRecalls(z)=recall;
        unsignedPrecisions(z)=precision;
        unsignedCountsMatrices{z}=countsMatrix;
    end
end