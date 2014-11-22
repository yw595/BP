function [fluxData fluxStds metsData metsStds] = readFluxAndMetsData(idxsToSkip1,idxsToSkip2)
    [junk1 junk2 totalData]=xlsread('../ChubukovData/inline-supplementary-material-22.xlsx',4);
    fluxData=[];
    for j=4:11
        for i=5:50
            fluxData(i-4,j-3)=totalData{i,j};
        end
    end
    fluxData(:,idxsToSkip1)=[];
    fluxData(idxsToSkip2,:)=[];
    
    %get standard deviations of 8 conditions besides glucose
    fluxStds=[];
    for j=13:20
        for i=5:50
            fluxStds(i-4,j-12)=totalData{i,j};
        end
    end
    fluxStds(:,idxsToSkip1)=[];
    fluxStds(idxsToSkip2,:)=[];
    
    [junk1 junk2 totalData]=xlsread('../ChubukovData/inline-supplementary-material-2.xlsx',3);
    
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
    
    %get standard deviations of 8 conditions besides glucose
    metsStds=[];
    for j=11:18
        for i=5:38
            metsStds(i-4,j-10)=totalData{i,j};
        end
    end
end