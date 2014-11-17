function averageModels(suffix,idxsToSkip1,idxsToSkip2)
    [metsNames fluxNames]=renameMetsAndRxns();
    fluxNamesIdxs=1:length(fluxNames);
    fluxNamesIdxs(idxsToSkip2)=[];
    fluxNames=fluxNames(fluxNamesIdxs);
    
    %take average value among all W matrices
    files=dir(suffix);
    WAverage=zeros(length(fluxNames),length(metsNames));
    numFilesRead=0;
    for i=1:length(files)
        if(sum(regexp(files(i).name,'\d.txt'))~=0)
            numFilesRead=numFilesRead+1;
            W=readMatrix([suffix '\' files(i).name]);
            WAverage=WAverage+W;
        end
    end
    WAverage=WAverage/numFilesRead;
    
    writeMatrix(WAverage,[suffix '\ModelAverage.txt'],0)
    writeMatrix(WAverage,[suffix '\ModelAverage.sif'],0,1,0,.25,metsNames,fluxNames)
    writeMatrix(WAverage,[suffix '\ModelAverageedges.txt'],0,0,1,.25,metsNames,fluxNames)
end