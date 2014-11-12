function averageModels(suffix,idxsToSkip1,idxsToSkip2)
    [metsNames fluxNames]=renameMetsAndRxns();
    fluxNamesIdxs=1:length(fluxNames);
    fluxNamesIdxs(idxsToSkip2)=[];
    fluxNames=fluxNames(fluxNamesIdxs);
    files=dir(suffix);
    WAverage=zeros(length(fluxNames),length(metsNames));
    numFilesRead=0;
    for i=1:length(files)
        if(sum(regexp(files(i).name,'\d.txt'))~=0)
            numFilesRead=numFilesRead+1;
            j=0;
            W=zeros(length(fluxNames),length(metsNames));
            modelFID=fopen([suffix '\' files(i).name]);
            line=fgetl(modelFID);
            j=j+1;
            while(line~=-1)
                words=strsplit(line,' ');
                words=words(2:end);
                W(j,1)=str2num(words{1});
                for k=2:length(words)
                    W(j,k)=str2num(words{k});
                end
                line=fgetl(modelFID);
                j=j+1;
            end
            fclose(modelFID);
            WAverage=WAverage+W;
        end
    end
    WAverage=WAverage/numFilesRead;
    ModelAverageTxtFID=fopen([suffix '\ModelAverage.txt'],'w');
    fprintf(ModelAverageTxtFID,'EdgeProbability\n');
    ModelAverageSifFID=fopen([suffix '\ModelAverage.sif'],'w');
    ModelAverageEdgesFID=fopen([suffix '\ModelAverageedges.txt'],'w');
    for i=1:size(WAverage,1)
        for j=1:size(WAverage,2)
            if(j==size(WAverage,2))
                fprintf(ModelAverageTxtFID,' %f\n',WAverage(i,j));
            else
                fprintf(ModelAverageTxtFID,' %f',WAverage(i,j));
            end
            if(abs(WAverage(i,j))>=.25)
                if(WAverage(i,j)<0)
                    fprintf(ModelAverageSifFID,'%s (inhibits) %s\n',metsNames{j},fluxNames{i});
                    fprintf(ModelAverageEdgesFID,'%s (inhibits) %s = %f\n',metsNames{j},fluxNames{i},WAverage(i,j));
                end
                if(WAverage(i,j)>0)
                    fprintf(ModelAverageSifFID,'%s (activates) %s\n',metsNames{j},fluxNames{i});
                    fprintf(ModelAverageEdgesFID,'%s (activates) %s = %f\n',metsNames{j},fluxNames{i},WAverage(i,j));
                end
            end
        end
    end
    fclose(ModelAverageTxtFID);
    fclose(ModelAverageSifFID);
    fclose(ModelAverageEdgesFID);
end