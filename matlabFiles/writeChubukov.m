function writeChubukov(suffix,idxsToSkip1, idxsToSkip2,localMin,globalMin,...
    perturbations,duplicate,influenceHardPrior,correlationHardPrior,...
    influenceSoftPrior,correlationSoftPrior,useFluxLog,useMetsLog,deciOption,stoichioSoftPrior,...
    useLinearAnsatz,beta,lambda)
    corrThresh=0.5;
    
    if(~exist(['../' suffix],'dir'))
        eval(['mkdir ../' suffix]);
    end
    
    [fluxData fluxStds metsData metsStds] = readFluxAndMetsData(idxsToSkip1,idxsToSkip2);

    [metsNames fluxNames]=renameMetsAndRxns();
    fluxNamesIdxs=1:length(fluxNames);
    fluxNamesIdxs(idxsToSkip2)=[];
    fluxNames=fluxNames(fluxNamesIdxs);

    fluxData = prepareData(fluxData,fluxStds,duplicate,localMin,globalMin,useFluxLog);
    writeMatrix(fluxData,['../' suffix '/fluxData.txt']);

    fluxPert = zeros(size(fluxData,1),size(fluxData,2));
    fluxPert(1:8,:)=atanh(fluxData(1:8,:)./(1.1*max(fluxData(1:8,:),[],2)*ones(1,size(fluxData,2))));
    writeMatrix(fluxPert,['../' suffix '/fluxPert.txt']);

    fluxNameFID=fopen(['../' suffix '/fluxName.txt'],'w');
    for i=1:size(fluxData,1)
        fprintf(fluxNameFID,'1 %s\n',fluxNames{i});
    end
    fclose(fluxNameFID);

    SMatrix=formSMatrix(metsNames,fluxNames);
    rxnInfluenceMatrix=formRxnInfluenceMatrix(metsNames,fluxNames,SMatrix);
    
    corrMatrix=makeCorrMatrix(fluxData,fluxData);
    sigCorrMatrix = corrMatrix;
    sigCorrMatrix(abs(corrMatrix)<corrThresh)=0;
    
    if(influenceSoftPrior)
        writeMatrix(reshapeMatrix(rxnInfluenceMatrix),['../' suffix '/fluxSoftPrior.txt'],1);
    elseif(correlationSoftPrior)
        writeMatrix(reshapeMatrix(sign(sigCorrMatrix)),['../' suffix '/fluxSoftPrior.txt'],1);
    end
    
    if(influenceHardPrior)
        fluxHardPrior=abs(rxnInfluenceMatrix);
    elseif(correlationHardPrior)
        fluxHardPrior=sigCorrMatrix;
    else
        fluxHardPrior=ones(length(fluxNames),length(fluxNames));
    end
    writeMatrix(fluxHardPrior,['../' suffix '/fluxHardPrior.txt'],1);

    metsData = prepareData(metsData,metsStds,duplicate,0,0,useMetsLog);    
    writeMatrix(metsData,['../' suffix '/metsData.txt']);
    
    if(stoichioSoftPrior)
        writeMatrix(reshapeMatrix(sign(SMatrix)),['../' suffix '/metsSoftPrior.txt'],1);
    end

    metsHardPrior=ones(size(SMatrix,1),size(SMatrix,2));
    if(stoichioSoftPrior)
        metsHardPrior=SMatrix~=0;
    end
    writeMatrix(metsHardPrior,['../' suffix '/metsHardPrior.txt'],1);
    
    metsNameFID=fopen(['../' suffix '/metsName.txt'],'w');
    for i=1:size(metsData,1)
        fprintf(metsNameFID,'%s\n',metsNames{i});
    end
    
    inputFID=fopen(['../' suffix '/input.txt'],'w');
    fprintf(inputFID,'%d\n',0);                         %doTestInd2Sub
    fprintf(inputFID,'%s\n',suffix);                    %inputDir
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
        fprintf(inputFID,'%d\n',sum(sum(sigCorrMatrix>0)));
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