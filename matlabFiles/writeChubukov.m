function writeChubukov(suffix,idxsToSkip1, idxsToSkip2,localMin,globalMin,...
    perturbations,duplicate,influenceHardPrior,correlationHardPrior,...
    influenceSoftPrior,correlationSoftPrior,useFluxLog,useMetsLog,deciOption,stoichioSoftPrior,...
    useLinearAnsatz,beta,lambda)
    corrThresh=0.5;
    
    [fluxData fluxStds metsData metsStds] = readFluxAndMetsData(idxsToSkip1,idxsToSkip2);

    [metsNames fluxNames]=renameMetsAndRxns();
    fluxNamesIdxs=1:length(fluxNames);
    fluxNamesIdxs(idxsToSkip2)=[];
    fluxNames=fluxNames(fluxNamesIdxs);

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
    
    corrMatrix=makeCorrMatrix(fluxData);

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

    metsData = prepareData(metsData,metsStds,duplicate,0,0,useMetsLog);    
    writeMatrix(metsData,['metsData' suffix '.txt']);
    
    if(stoichioSoftPrior)
        writeMatrix(reshapeMatrix(sign(SMatrix)),['metsSoftPrior' suffix '.txt']);
    end

    metsHardPrior=ones(size(SMatrix,1),size(SMatrix,2));
    if(stoichioSoftPrior)
        metsHardPrior=SMatrix~=0;
    end
    writeMatrix(metsHardPrior,['metsHardPrior' suffix '.txt'],1);
    
    metsNameFID=fopen(['metsName' suffix '.txt'],'w');
    for i=1:size(metsData,1)
        fprintf(metsNameFID,'%s\n',metsNames{i});
    end
    
    inputFID=fopen(['input' suffix '.txt'],'w');
    fprintf(inputFID,'%s\n',suffix);                    %Session ID
    if(duplicate)                                       %Nexpts
        fprintf(inputFID,'%d\n',size(fluxData,2)/5);
    else
        fprintf(inputFID,'%d\n',size(fluxData,2));
    end
    fprintf(inputFID,'%d\n',size(fluxData,1));          %Nnodes
    fprintf(inputFID,'%d\n',15);                        %Nwvals
    fprintf(inputFID,'%s\n','2.00');                    %MaxWvalss
    fprintf(inputFID,'%s\n','1.00E-06');                %thresh
    fprintf(inputFID,'%s\n',num2str(lambda));           %lambda
    fprintf(inputFID,'%s\n',num2str(beta));             %beta
    if(influenceSoftPrior)                              %Nprior
        fprintf(inputFID,'%d\n',sum(sum(rxnInfluenceMatrix~=0)));
    elseif(correlationSoftPrior)
        fprintf(inputFID,'%d\n',sum(sum(abs(corrMatrix)>=corrThresh)));
    elseif(stoichioSoftPrior)
        fprintf(inputFID,'%d\n',sum(sum(SMatrix~=0)));
    else
        fprintf(inputFID,'%d\n',0);
    end
    if(perturbations)                                   %Nobs
        fprintf(inputFID,'%d\n',size(fluxData,1)-8);
    else
        fprintf(inputFID,'%d\n',size(fluxData,1));
    end
    if(deciOption==1)                                   %Ndeci
        fprintf(inputFID,'%d\n',100);
    end
    if(deciOption==2 || deciOption==3)
        fprintf(inputFID,'%d\n',1);
    end
    fprintf(inputFID,'%d\n',10);                        %opt_maxSkip
    fprintf(inputFID,'%d\n',size(metsData,1));          %Nmets
    fprintf(inputFID,'%d\n',deciOption);                %deciOption
    fprintf(inputFID,'%d\n',size(metsData,1));          %Nparams
    fprintf(inputFID,'%d\n',useLinearAnsatz);           %useLinearAnsatz
    fclose(inputFID);
end