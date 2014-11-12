function readMarginals(suffix,modelFile,idxsToSkip2)
    [metsNames fluxNames]=renameMetsAndRxns();
    Nwvals=15;
    fluxNamesIdxs=1:length(fluxNames);
    fluxNamesIdxs(idxsToSkip2)=[];
    fluxNames=fluxNames(fluxNamesIdxs);
    interactionsMatrix=zeros(length(fluxNames),length(metsNames),Nwvals);
    for i=1:length(fluxNames)
        target=fluxNames{i};
        targetFID=fopen([suffix '\Marginals\' num2str(i) '_1_' target '.txt']);
        line=fgetl(targetFID);
        lineNum=1;
        while(line~=-1)
            words=strsplit(line,'  ');
            words=words(2:end);
            for j=1:length(words)
                interactionsMatrix(i,lineNum,j)=str2num(words{j});
            end
            line=fgetl(targetFID);
            lineNum=lineNum+1;
        end
        fclose(targetFID);
    end
    
    modelFID=fopen(modelFile);
    line=fgetl(modelFID);
    marginals=[];
    labels={};
    while(line~=-1)
        words=strsplit(line,' ');
        met=words{1};flux=words{3};
        marginals(end+1,:)=squeeze(interactionsMatrix(strcmp(fluxNames,flux),strcmp(metsNames,met),:))';
        labels{end+1}=[met ' ' flux];
        line=fgetl(modelFID);
    end
    fclose(modelFID);
    
    hIndex=1;
    while(hIndex*40<=size(marginals,1))
        if((hIndex-1)*40+40<=size(marginals,1))
            h=HeatMap(marginals((hIndex-1)*40+1:(hIndex-1)*40+40,:),...
            'RowLabels',labels((hIndex-1)*40+1:(hIndex-1)*40+40));
        else
            h=HeatMap(marginals((hIndex-1)*40+1:end,:),...
            'RowLabels',labels((hIndex-1)*40+1:end));
        end
        saveas(plot(h),[suffix '\Marginals' num2str(hIndex) '.png'])
        hIndex=hIndex+1;
    end
    close all force
end