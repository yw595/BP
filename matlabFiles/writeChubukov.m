function writeChubukov(suffix,idxsToSkip1, idxsToSkip2,localMin,globalMin,...
    perturbations,duplicate,influenceHardPrior,correlationHardPrior,...
    influenceSoftPrior,correlationSoftPrior,useFluxLog,useMetsLog,deciOption,stoichioSoftPrior,...
    useLinearAnsatz,beta,lambda)
    corrThresh=0.5;
    [junk1 junk2 totalData]=xlsread('../Chubukov/inline-supplementary-material-22.xlsx',4);

    %get reaction fluxNames in storefluxNames, take out spaces for fluxNames
    [metsNames fluxNames]=renameMetsAndRxns();
    fluxNamesIdxs=1:length(fluxNames);
    fluxNamesIdxs(idxsToSkip2)=[];
    fluxNames=fluxNames(fluxNamesIdxs);

    %get fluxData as all fluxes, fluxDataNorm as first condition
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

    fluxData = prepareData(fluxData,fluxStds,duplicate,localMin,globalMin,useFluxLog);
    writeMatrix(fluxData,['fluxData' suffix '.txt']);

    fluxPert = zeros(size(fluxData,1),size(fluxData,2));
    fluxPert(1:8,:)=fluxData(1:8,:)./(1.1*max(fluxData(1:8,:),2)*ones(8,size(fluxData,2)));
    writeMatrix(fluxPert,['fluxPert' suffix '.txt']);

    fluxNameFID=fopen(['fluxName' suffix '.txt'],'w');
    for i=1:size(fluxData,1)
        fprintf(fluxNameFID,'1 %s\n',fluxNames{i});
    end
    fclose(fluxNameFID);

    SMatrix=formSMatrix(metsNames,fluxNames);
    rxnInfluenceMatrix=formInfluenceMatrix(fluxNames,SMatrix,metsNames);
    
    corrMatrix=zeros(length(fluxNames),length(fluxNames));
    for i=1:length(fluxNames)
        for j=1:length(fluxNames)
            if(i~=j)
                row1=fluxData(i, fluxData(i,:) & fluxData(j,:) );
                row2=fluxData(j, fluxData(i,:) & fluxData(j,:) );
                if(length(row1)>3)
                    tempCorr=corr([row1' row2'],'type','Spearman');
                    corrMatrix(i,j)=tempCorr(1,2);
                end
            end
        end
    end

    if(influenceSoftPrior)
        writeMatrix(reshapeMatrix(rxnInfluenceMatrix),['fluxSoftPrior' suffix '.txt']);
    elseif(correlationSoftPrior)
        writeMatrix(reshapeMatrix(sign(corrMatrix (abs(corrMatrix)>=corrThresh) )),['fluxSoftPrior' suffix '.txt']);
    end
    
    if(influenceHardPrior)
        fluxHardPrior=abs(rxnInfluenceMatrix);
    elseif(correlationHardPrior)
        fluxHardPrior=abs(corrMatrix)>=corrThresh;
    else
        fluxHardPrior=ones(length(fluxNames),length(fluxNames));
    end
    writeMatrix(fluxHardPrior,['fluxHardPrior' suffix '.txt'],1);
    
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
    
    %get standard deviations of 8 conditions besides glucose
    metsStds=[];
    for j=11:18
        for i=5:38
            metsStds(i-4,j-10)=totalData{i,j};
        end
    end

    metsData = prepareData(metsData,metsStds,duplicate,0,0,useMetsLog);    
    writeMatrix(metsData,['metsData' suffix '.txt']);
    
    if(stoichioSoftPrior)
        writeMatrix(reshapeMatrix(sign(SMatrix)),['metsSoftPrior' suffix '.txt']);
    end

    metsHardPrior=ones(size(SMatrix,1),size(SMatrix,2));
    if(stoichioSoftPrior)
        metsHardPrior=SMatrix(i,j)~=0;
    end
    writeMatrix(metsHardPrior,['metsHardPrior' suffix '.txt'],1);
    
    metsNameFID=fopen(['metsName' suffix '.txt'],'w');
    for i=1:size(metsData,1)
        fprintf(metsNameFID,'%s\n',metsNames{i});
    end
    
    inputFID=fopen(['input' suffix '.txt'],'w');
    %Session ID
    fprintf(inputFID,'%s\n',suffix);
    %Nexpts
    if(duplicate)
        fprintf(inputFID,'%d\n',size(fluxData,2)/5);
    else
        fprintf(inputFID,'%d\n',size(fluxData,2));
    end
    %Nnodes
    fprintf(inputFID,'%d\n',size(fluxData,1));
    %Nwvals
    fprintf(inputFID,'%d\n',15);
    %MaxWvals
    fprintf(inputFID,'%s\n','2.00');
    %thresh
    fprintf(inputFID,'%s\n','1.00E-06');
    %lambda
    fprintf(inputFID,'%s\n',num2str(lambda));
    %beta
    fprintf(inputFID,'%s\n',num2str(beta));
    %Nprior
    if(influenceSoftPrior)
        fprintf(inputFID,'%d\n',sum(sum(rxnInfluenceMatrix~=0)));
    elseif(correlationSoftPrior)
        fprintf(inputFID,'%d\n',sum(sum(abs(corrMatrix)>=corrThresh)));
    elseif(stoichioSoftPrior)
        fprintf(inputFID,'%d\n',sum(sum(SMatrix~=0)));
    else
        fprintf(inputFID,'%d\n',0);
    end
    %Nobs
    if(perturbations)
        fprintf(inputFID,'%d\n',size(fluxData,1)-8);
    else
        fprintf(inputFID,'%d\n',size(fluxData,1));
    end
    %Ndeci
    if(deciOption==1)
        fprintf(inputFID,'%d\n',100);
    end
    if(deciOption==2 || deciOption==3)
        fprintf(inputFID,'%d\n',1);
    end
    %opt_maxSkip
    fprintf(inputFID,'%d\n',10);
    %Nmets
    fprintf(inputFID,'%d\n',size(metsData,1));
    %deciOption
    fprintf(inputFID,'%d\n',deciOption);
    %Nparams
    fprintf(inputFID,'%d\n',size(metsData,1));
    %useLinearAnsatz
    fprintf(inputFID,'%d\n',useLinearAnsatz);
    fclose(inputFID);
end