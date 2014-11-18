function [signedRecall signedPrecision signedCounts unsignedRecall unsignedPrecision unsignedCounts]= ...
    validateCorrelate(dataFile,idxsToSkip1,idxsToSkip2,isMets)
    corrThresh=0.5;
    
    [junk1 junk2 totalData]=xlsread('../Chubukov/inline-supplementary-material-22.xlsx',4);

    %get reaction fluxNames in storefluxNames, take out spaces for fluxNames
    [metsNames fluxNames]=renameMetsAndRxns();
    fluxNamesIdxs=1:length(fluxNames);
    fluxNamesIdxs(idxsToSkip2)=[];
    fluxNames=fluxNames(fluxNamesIdxs);

    %get fluxData as all fluxes, fluxDataNorm as first condition
    fluxDataNorm=[];
    fluxData=[];
    for j=4:11
        for i=5:50
            fluxData(i-4,j-3)=totalData{i,j};
        end
    end
    fluxData(:,idxsToSkip1)=[];
    fluxData(idxsToSkip2,:)=[];

    fluxCorrMatrix=zeros(length(fluxNames),length(fluxNames));
    fluxCorrs=[];
    for i=1:length(fluxNames)
        for j=1:length(fluxNames)
            if(i~=j)
                row1=[];
                row2=[];
                for k=1:size(fluxData,2)
                    if(abs(fluxData(i,k))>0 && abs(fluxData(j,k))>0)
                        row1(end+1)=fluxData(i,k);
                        row2(end+1)=fluxData(j,k);
                    end
                end
                if(length(row1)>3)
                    tempCorr=corr([row1' row2'],'type','Spearman');
                    fluxCorrMatrix(i,j)=tempCorr(1,2);
                    if(i>j)
                        fluxCorrs(end+1,1)=tempCorr(1,2);
                        fluxCorrs(end,2)=i;
                        fluxCorrs(end,3)=j;
                    end
                end
            end
        end
    end
    fluxCorrsTemp=[];
    for i=1:length(fluxCorrs)
        if(abs(fluxCorrs(i,1))>0.7)
            fluxCorrsTemp(end+1,:)=fluxCorrs(i,:);
        end
    end
    fluxCorrs=fluxCorrsTemp;
    length(fluxCorrs)
    %[junk sortIdxs]=sort(abs(fluxCorrs(:,1)));
    %fluxCorrs=fluxCorrs(sortIdxs,:);
    fluxCorrMatrix=zeros(size(fluxCorrMatrix,1),size(fluxCorrMatrix,2));
    for i=1:length(fluxCorrs)
        fluxCorrMatrix(fluxCorrs(i,2),fluxCorrs(i,3))=fluxCorrs(i,1);
    end
    %sum(sum(fluxCorrMatrix))
    
    [junk1 junk2 totalData]=xlsread('../Chubukov/inline-supplementary-material-2.xlsx',3);
    
    metsData=[];
    for j=2:9
        for i=5:38
            if(strcmp(totalData{i,j},'n/a'))
                metsData(i-4,j-1)=0;
            else
                metsData(i-4,j-1)=totalData{i,j};
            end
        end
    end
    metsData(:,idxsToSkip1)=[];
    
    fluxMetsCorrMatrix=zeros(length(metsNames),length(fluxNames));
    fluxMetsCorrs=[];
    for i=1:length(metsNames)
        for j=1:length(fluxNames)
            row1=[];
            row2=[];
            for k=1:size(fluxData,2)
                if(abs(fluxData(j,k))>0 && abs(metsData(i,k))>0)
                    row1(end+1)=fluxData(j,k);
                    row2(end+1)=metsData(i,k);
                end
            end
            if(length(row1)>3)
                tempCorr=corr([row1' row2'],'type','Spearman');
                fluxMetsCorrMatrix(i,j)=tempCorr(1,2);
                fluxMetsCorrs(end+1,1)=tempCorr(1,2);
                fluxMetsCorrs(end,2)=i;
                fluxMetsCorrs(end,3)=j;
            end
        end
    end
    
    [junk sortIdxs]=sort(abs(fluxMetsCorrs(:,1)));
    fluxMetsCorrs=fluxMetsCorrs(sortIdxs(end-99:end),:);
    fluxMetsCorrMatrix=zeros(size(fluxMetsCorrMatrix,1),size(fluxMetsCorrMatrix,2));
    for i=1:length(fluxMetsCorrs)
        fluxMetsCorrMatrix(fluxMetsCorrs(i,2),fluxMetsCorrs(i,3))=fluxMetsCorrs(i,1);
    end
    
    nodeInfluenceMatrix=[];
    nodeNames={};
    if(isMets)
        %fluxMetsCorrMatrix(abs(fluxMetsCorrMatrix)<=corrThresh)=0;
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

    %from edges file, read in sign of influence in InfluenceMatrix, size of
    %influence in InfluenceMatrixContinuous
    WNodeInfluenceMatrix=zeros(size(nodeInfluenceMatrix,1),size(nodeInfluenceMatrix,2));
    WNodeInfluenceMatrixContinuous=zeros(size(nodeInfluenceMatrix,1),size(nodeInfluenceMatrix,2));
    WNodeEdges=[];
    WFID=fopen(dataFile);
    line=fgetl(WFID);
    while(line~=-1)
        firstNode=line(1:regexp(line,'(')-2);
        secondNode=line(regexp(line,')')+2:regexp(line,'=')-2);
        influence=line(regexp(line,'(')+1:regexp(line,')')-1);
        if(sum(strcmp(secondNode,nodeNames))~=0 && sum(strcmp(firstNode,nodeNames))~=0)
            WNodeInfluenceMatrixContinuous(strcmp(nodeNames,firstNode),strcmp(nodeNames,secondNode))=...
            str2num(line(min(regexp(line,'[-|\d]+\.\d+')):end));
            WNodeEdges(end+1,1)=sign(str2num(line(min(regexp(line,'[-|\d]+\.\d+')):end)));
            WNodeEdges(end,2)=find(strcmp(firstNode,nodeNames));
            WNodeEdges(end,3)=find(strcmp(secondNode,nodeNames));
            if(strcmp(influence,'activates'))
                WNodeInfluenceMatrix(strcmp(nodeNames,firstNode),strcmp(nodeNames,secondNode))=1;
            elseif(strcmp(influence,'inhibits'))
                WNodeInfluenceMatrix(strcmp(nodeNames,firstNode),strcmp(nodeNames,secondNode))=-1;
            end
        end
        line=fgetl(WFID);
    end
    %WNodeEdgesTemp=[];
    %for i=1:length(WNodeEdges)
        %if(abs(WNodeEdges(i,1))>0.8)
            %WNodeEdgesTemp(end+1,:)=WNodeEdges(i,:);
        %end
    %end
    %WNodeEdges=WNodeEdgesTemp;
    [junk sortIdxs]=sort(abs(WNodeEdges(:,1)));
    WNodeEdges=WNodeEdges(sortIdxs,:);
    WNodeInfluenceMatrix=zeros(size(WNodeInfluenceMatrix,1),size(WNodeInfluenceMatrix,2));
    for i=1:length(WNodeEdges)
        WNodeInfluenceMatrix(WNodeEdges(i,2),WNodeEdges(i,3))=WNodeEdges(i,1);
    end

    diffMatrixFID=fopen('diffMatrix.sif','w');
    for i=1:size(WNodeInfluenceMatrix,1)
        for j=1:size(WNodeInfluenceMatrix,2)
            if(WNodeInfluenceMatrix(i,j)~=nodeInfluenceMatrix(i,j) && WNodeInfluenceMatrix(i,j)>0)
                fprintf(diffMatrixFID,'%s activates %s\n',nodeNames{i},nodeNames{j});
            elseif(WNodeInfluenceMatrix(i,j)~=nodeInfluenceMatrix(i,j) && WNodeInfluenceMatrix(i,j)<0)
                fprintf(diffMatrixFID,'%s inhibitss %s\n',nodeNames{i},nodeNames{j});
            end
        end
    end
    fclose(diffMatrixFID);
    
    diffMatrixEdgesFID=fopen('diffMatrixEdges.txt','w');
    fprintf(diffMatrixEdgesFID,'EdgeProbability\n');
    for i=1:size(WNodeInfluenceMatrix,1)
        for j=1:size(WNodeInfluenceMatrix,2)
            if(WNodeInfluenceMatrix(i,j)~=nodeInfluenceMatrix(i,j) && WNodeInfluenceMatrix(i,j)>0)
                fprintf(diffMatrixEdgesFID,'%s (activates) %s = %f\n',nodeNames{i},nodeNames{j},WNodeInfluenceMatrix(i,j));
            elseif(WNodeInfluenceMatrix(i,j)~=nodeInfluenceMatrix(i,j) && WNodeInfluenceMatrix(i,j)<0)
                fprintf(diffMatrixEdgesFID,'%s (inhibits) %s = %f\n',nodeNames{i},nodeNames{j},WNodeInfluenceMatrix(i,j));
            end
        end
    end
    fclose(diffMatrixEdgesFID);
    
    signedRecall=[];signedPrecision=[];unsignedRecall=[];unsignedPrecision=[];signedCounts={};unsignedCounts={};
    for z=1:10 
        %form NodeInfluenceMatrix with all influences up to distance of cutoff
        %allowed
        cutoff=z;
        for i=1:length(nodeNames)
            for j=1:length(nodeNames)
                [junk1,path,junk2]=graphshortestpath(sparse(absNodeInfluenceMatrix),i,j);
                if(length(path)<=cutoff+1 && length(path)>1)
                    pathSign=1;
                    for k=1:length(path)-1
                        pathSign=pathSign*storeNodeInfluenceMatrix(path(k),path(k+1));
                    end
                    nodeInfluenceMatrix(i,j)=pathSign;
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
        for i=1:length(nodeNames)
            for j=1:length(nodeNames)
                ijthPossibility=possibilities{storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2};
                countsMatrix(storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)= ...
                countsMatrix(storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)+1;
                %if false positive or negative, use extended NodeInfluence instead,
                %adjusting countsMatrix. If still false, write out Nodes, weight of
                %predicted interaction, and correct interaction
                if(strcmp(ijthPossibility,'False Positive') || strcmp(ijthPossibility,'False Negative'))
                    ijthPossibility=possibilities{nodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2};
                    countsMatrix(storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)= ...
                    countsMatrix(storeNodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)-1;
                    countsMatrix(nodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)= ...
                    countsMatrix(nodeInfluenceMatrix(i,j)+2,WNodeInfluenceMatrix(i,j)+2)+1;
                end
            end
        end
    
        recall=(countsMatrix(1,1)+countsMatrix(3,3))/(sum(countsMatrix(1,:)+countsMatrix(3,:)));
        precision=(countsMatrix(1,1)+countsMatrix(3,3))/(sum(countsMatrix(:,1)+countsMatrix(:,3)));
        signedRecall(z)=recall;
        signedPrecision(z)=precision;
        signedCounts{z}=countsMatrix;
        
        %iterate over NodeInfluenceMatrix and W, calculate accuracy statistics
        possibilities=cell(2,2);
        possibilities{1,1}='True Positive';
        possibilities{1,2}='False Negative';
        possibilities{2,1}='False Positive';
        possibilities{2,2}='True Negative';
        countsMatrix=zeros(2,2);
        for i=1:length(nodeNames)
            for j=1:length(nodeNames)
                storeCoordinate=abs(storeNodeInfluenceMatrix(i,j));
                if(storeCoordinate==0)
                    storeCoordinate=2;
                end
                WCoordinate=abs(WNodeInfluenceMatrix(i,j));
                if(WCoordinate==0)
                    WCoordinate=2;
                end

                ijthPossibility=possibilities{storeCoordinate,WCoordinate};
                countsMatrix(storeCoordinate,WCoordinate)= ...
                countsMatrix(storeCoordinate,WCoordinate)+1;
                %if false positive or negative, use extended NodeInfluence instead,
                %adjusting countsMatrix. If still false, write out Nodes, weight of
                %predicted interaction, and correct interaction
                if(strcmp(ijthPossibility,'False Positive') || strcmp(ijthPossibility,'False Negative'))
                    ijthPossibility=possibilities{storeCoordinate,WCoordinate};
                    countsMatrix(storeCoordinate,WCoordinate)= ...
                    countsMatrix(storeCoordinate,WCoordinate)-1;
                    
                    storeCoordinate=abs(nodeInfluenceMatrix(i,j));
                    if(storeCoordinate==0)
                        storeCoordinate=2;
                    end
                    WCoordinate=abs(WNodeInfluenceMatrix(i,j));
                    if(WCoordinate==0)
                        WCoordinate=2;
                    end
                    
                    countsMatrix(storeCoordinate,WCoordinate)= ...
                    countsMatrix(storeCoordinate,WCoordinate)+1;
                end
            end
        end
    
        recall=countsMatrix(1,1)/sum(countsMatrix(1,:));
        precision=countsMatrix(1,1)/sum(countsMatrix(:,1));
        unsignedRecall(z)=recall;
        unsignedPrecision(z)=precision;
        unsignedCounts{z}=countsMatrix;
    end
end