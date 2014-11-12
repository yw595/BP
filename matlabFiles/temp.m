corrThreshold=0.9;

[junk1 junk2 totalData]=xlsread('../Chubukov/inline-supplementary-material-22.xlsx',4);

%get reaction fluxNames in storefluxNames, take out spaces for fluxNames
[metsNames fluxNames]=renameMetsAndRxns();
%fluxNamesIdxs=1:length(fluxNames);
%fluxNamesIdxs(idxsToSkip2)=[];
%fluxNames=fluxNames(fluxNamesIdxs);

%get fluxData as all fluxes, fluxDataNorm as first condition
fluxDataNorm=[];
fluxData=[];
for j=4:11
    for i=5:50
        fluxData(i-4,j-3)=totalData{i,j};
    end
end

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
                %if(length(fluxCorrs)<100)
                if(i>j)
                    fluxCorrs(end+1,1)=tempCorr(1,2);
                    fluxCorrs(end,2)=i;
                    fluxCorrs(end,3)=j;
                end
                %end
            end
        end
    end
end

[junk sortIdxs]=sort(abs(fluxCorrs(:,1)));
fluxCorrs=fluxCorrs(sortIdxs(end-99:end),:);

%write out storeRxnInfluence as .sif file
if(0)
corrNetworkSif=fopen('corrNetwork.sif','w');
for i=1:length(fluxCorrs)
    %for j=1:length(fluxNames)
        %if(fluxCorrMatrix(i,j)<=-corrThreshold)
            fprintf(corrNetworkSif,[fluxNames{i} ' inhibits ' fluxNames{j} '\n']);
        %elseif(fluxCorrMatrix(i,j)>=corrThreshold)
            fprintf(corrNetworkSif,[fluxNames{i} ' activates ' fluxNames{j} '\n']);
        %end
    %end
end
fclose(corrNetworkSif);
end

corrNetworkSif=fopen('corrNetwork.sif','w');
for i=1:length(fluxCorrs)
    %for j=1:length(fluxNames)
        if(fluxCorrs(i)<0)
            fprintf(corrNetworkSif,[fluxNames{fluxCorrs(i,2)} ' inhibits ' fluxNames{fluxCorrs(i,3)} '\n']);
        elseif(fluxCorrs(i)>0)
            fprintf(corrNetworkSif,[fluxNames{fluxCorrs(i,2)} ' activates ' fluxNames{fluxCorrs(i,3)} '\n']);
        end
    %end
end
fclose(corrNetworkSif);

corrNetworkEdges=fopen('corrNetworkEdges.txt','w');
fprintf(corrNetworkEdges,'EdgeProbability\n');
for i=1:length(fluxCorrs)
    %for j=1:length(fluxNames)
        if(fluxCorrs(i)<0)
            fprintf(corrNetworkEdges,[fluxNames{fluxCorrs(i,2)} ' (inhibits) ' fluxNames{fluxCorrs(i,3)} ' = ' num2str(fluxCorrs(i,1)) '\n']);
        elseif(fluxCorrs(i)>0)
            fprintf(corrNetworkEdges,[fluxNames{fluxCorrs(i,2)} ' (activates) ' fluxNames{fluxCorrs(i,3)} ' = ' num2str(fluxCorrs(i,1)) '\n']);
        end
    %end
end
fclose(corrNetworkEdges);

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
%metsData(:,idxsToSkip1)=[];

corrThreshold=0.6;

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
            %if(length(fluxMetsCorrs)<100)
            %if(i>j)
                fluxMetsCorrs(end+1,1)=tempCorr(1,2);
                fluxMetsCorrs(end,2)=i;
                fluxMetsCorrs(end,3)=j;
            %end
            %end
        end
    end
end

[junk sortIdxs]=sort(abs(fluxMetsCorrs(:,1)));
fluxMetsCorrs=fluxMetsCorrs(sortIdxs(end-99:end),:);

%write out storeRxnInfluence as .sif file
if(0)
corrMetsNetworkSif=fopen('corrMetsNetwork.sif','w');
for j=1:length(fluxNames)
    for i=1:length(metsNames)
        if(fluxMetsCorrMatrix(i,j)<=-corrThreshold)
            fprintf(corrMetsNetworkSif,[fluxNames{j} ' inhibits ' metsNames{i} '\n']);
        elseif(fluxMetsCorrMatrix(i,j)>=corrThreshold)
            fprintf(corrMetsNetworkSif,[fluxNames{j} ' activates ' metsNames{i} '\n']);
        end
    end
end
fclose(corrMetsNetworkSif);
end

corrMetsNetworkSif=fopen('corrMetsNetwork.sif','w');
for i=1:length(fluxMetsCorrs)
    %for j=1:length(fluxNames)
        if(fluxMetsCorrs(i)<0)
            fprintf(corrMetsNetworkSif,[metsNames{fluxMetsCorrs(i,2)} ' inhibits ' fluxNames{fluxMetsCorrs(i,3)} '\n']);
        elseif(fluxMetsCorrs(i)>0)
            fprintf(corrMetsNetworkSif,[metsNames{fluxMetsCorrs(i,2)} ' activates ' fluxNames{fluxMetsCorrs(i,3)} '\n']);
        end
    %end
end
fclose(corrMetsNetworkSif);

corrMetsNetworkEdges=fopen('corrMetsNetworkEdges.txt','w');
fprintf(corrMetsNetworkEdges,'EdgeProbability\n');
for i=1:length(fluxMetsCorrs)
    %for j=1:length(fluxNames)
        if(fluxMetsCorrs(i)<0)
            fprintf(corrMetsNetworkEdges,[metsNames{fluxMetsCorrs(i,2)} ' (inhibits) ' fluxNames{fluxMetsCorrs(i,3)} ' = ' num2str(fluxMetsCorrs(i,1)) '\n']);
        elseif(fluxMetsCorrs(i)>0)
            fprintf(corrMetsNetworkEdges,[metsNames{fluxMetsCorrs(i,2)} ' (activates) ' fluxNames{fluxMetsCorrs(i,3)} ' = ' num2str(fluxMetsCorrs(i,1)) '\n']);
        end
    %end
end
fclose(corrMetsNetworkEdges);