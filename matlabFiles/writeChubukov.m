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
    fluxDataNorm=[];
    fluxData=[];
    for j=4:11
        for i=5:50
            fluxData(i-4,j-3)=totalData{i,j};
        end
    end
    fluxData(:,idxsToSkip1)=[];
    fluxData(idxsToSkip2,:)=[];
    
    %get standard deviations of 8 conditions besides glucose
    fluxStdevs=[];
    for j=13:20
        for i=5:50
            fluxStdevs(i-4,j-12)=totalData{i,j};
        end
    end
    fluxStdevs(:,idxsToSkip1)=[];
    fluxStdevs(idxsToSkip2,:)=[];

    %use standard deviations to replicate fluxData 5 times with Gaussian noise
    fluxDataTemp=[];
    for i=1:size(fluxData,1)
        for j=1:size(fluxData,2)
            for k=1:5
                fluxDataTemp(i,(j-1)*5+k)=fluxData(i,j)+fluxStdevs(i,j)*randn(1,1);
            end
        end
    end
    if(duplicate)
        fluxData=fluxDataTemp;
    end

    %calculate minimum absolute flux, subtract from all fluxes so all fluxes
    %are positive, thus fluxDataNorm values will be somewhere in that positive
    %range
    if(localMin)
        for i=1:size(fluxData,1)
            minFluxDataVal=min(fluxData(i,:));
            fluxData(i,:)=fluxData(i,:)-minFluxDataVal;
        end
    end
    if(globalMin)
        minFluxDataVal=min(min(fluxData));
        fluxData=fluxData-minFluxDataVal;
    end
    epsFlux=min(min(abs(fluxData(fluxData~=0))));
    if(useFluxLog)
        fluxDataNorm=fluxData(:,1);
        if(duplicate)
            fluxData=fluxData(:,6:end);
        else
            fluxData=fluxData(:,2:end);
        end
        for i=1:size(fluxData,2)
            for j=1:size(fluxData,1)
                if(fluxData(j,i)==0 || fluxDataNorm(j)==0)
                    fluxData(j,i)=log((fluxData(j,i)+epsFlux)/(fluxDataNorm(j)+epsFlux));
                else
                    fluxData(j,i)=log(fluxData(j,i)/fluxDataNorm(j));
                end
            end
        end
    end

    %write out log-fold values, setting limits at abs(99) so number of digits
    %to left of decimal <=2
    fluxDataFID=fopen(['fluxData' suffix '.txt'],'w');
    for i=1:size(fluxData,1)
        for j=1:size(fluxData,2)
            if(j==size(fluxData,2))
                fprintf(fluxDataFID,'%2.2f\n',fluxData(i,j));
            else
                fprintf(fluxDataFID,'%2.2f ',fluxData(i,j));
            end
        end
    end
    fclose(fluxDataFID);

    fluxPertFID=fopen(['fluxPert' suffix '.txt'],'w');
    for i=1:size(fluxData,1)
        for j=1:size(fluxData,2)
            writeValue=0;
            if(perturbations && i<=8 && fluxData(i,j)~=0)
                writeValue=fluxData(i,j)/(1.1*max(fluxData(i,:)));
            end
            if(j==size(fluxData,2))
                fprintf(fluxPertFID,'%2.2f\n',writeValue);
            else
                fprintf(fluxPertFID,'%2.2f ',writeValue);
            end
        end
    end
    fclose(fluxPertFID);

    %write out fluxNames without spaces so BP can read them
    fluxNameFID=fopen(['fluxName' suffix '.txt'],'w');
    for i=1:size(fluxData,1)
        fprintf(fluxNameFID,'1 %s\n',fluxNames{i});
    end
    fclose(fluxNameFID);

    %parse equations written out in rxnEqs.txt
    SMatrix=formSMatrix(metsNames,fluxNames);

    %form rxnInfluenceMatrix by dot-product of S-matrix columns for reaction,
    %filter out elements at CO2 rows, and remove Fba and Tpi interactions since
    %are ambiguous
    rxnInfluenceMatrix=formInfluenceMatrix(fluxNames,SMatrix,metsNames);
    
    corrMatrix=zeros(length(fluxNames),length(fluxNames));
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
                    corrMatrix(i,j)=tempCorr(1,2);
                end
            end
        end
    end

    %write out all entries in rxnInfluenceMatrix as prior
    fluxSoftPriorFID=fopen(['fluxSoftPrior' suffix '.txt'],'w');
    for i=1:length(fluxNames)
        for j=1:length(fluxNames)
            if(influenceSoftPrior)
                if(rxnInfluenceMatrix(i,j)~=0)
                    fprintf(fluxSoftPriorFID,'%d %d %d\n',i,j,rxnInfluenceMatrix(i,j));
                end
            elseif(correlationSoftPrior)
                if(abs(corrMatrix(i,j))>=corrThresh)
                    fprintf(fluxSoftPriorFID,'%d %d %d\n',i,j,sign(corrMatrix(i,j)));
                end
            end
        end
    end
    fclose(fluxSoftPriorFID);
    
    fluxHardPriorFID=fopen(['fluxHardPrior' suffix '.txt'],'w');
    for i=1:length(fluxNames)
        for j=1:length(fluxNames)
            writeValue=0;
            if(influenceHardPrior)
                writeValue=abs(rxnInfluenceMatrix(i,j));
            elseif(correlationHardPrior)
                writeValue=abs(corrMatrix(i,j))>=corrThresh;
            else
                writeValue=1;
            end
            if(j==length(fluxNames))
                fprintf(fluxHardPriorFID,'%d\n',writeValue);
            else
                fprintf(fluxHardPriorFID,'%d ',writeValue);
            end
        end
    end
    fclose(fluxHardPriorFID);
    
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
    metsStdevs=[];
    for j=11:18
        for i=5:38
            metsStdevs(i-4,j-10)=totalData{i,j};
        end
    end

    %use standard deviations to replicate fluxData 5 times with Gaussian noise
    metsDataTemp=[];
    for i=1:size(metsData,1)
        for j=1:size(metsData,2)
            for k=1:5
                metsDataTemp(i,(j-1)*5+k)=metsData(i,j)+metsStdevs(i,j)*randn(1,1);
            end
        end
    end
    if(duplicate)
        metsData=metsDataTemp;
    end
    
    epsMets=min(min(abs(metsData(metsData~=0))));
    
    if(useMetsLog)
        metsDataNorm=metsData(:,1);
        if(duplicate)
            metsData=metsData(:,6:end);
        else
            metsData=metsData(:,2:end);
        end
        for i=1:size(metsData,2)
            for j=1:size(metsData,1)
                if(metsData(j,i)==0 || metsDataNorm(j)==0)
                    metsData(j,i)=log((metsData(j,i)+epsFlux)./(metsDataNorm(j)+epsFlux));
                else
                    metsData(j,i)=log(metsData(j,i)/metsDataNorm(j));
                end
            end
        end
    end
    
    metsDataFID=fopen(['metsData' suffix '.txt'],'w');
    for i=1:size(metsData,1)
        for j=1:size(metsData,2)
            if(j==size(metsData,2))
                fprintf(metsDataFID,'%2.2f\n',metsData(i,j));
            else
                fprintf(metsDataFID,'%2.2f ',metsData(i,j));
            end
        end
    end
    fclose(metsDataFID);
    
    metsSoftPriorFID=fopen(['metsSoftPrior' suffix '.txt'],'w');
    for i=1:size(SMatrix,1)
        for j=1:size(SMatrix,2)
            if(stoichioSoftPrior)
                if(SMatrix(i,j)~=0)
                    fprintf(metsSoftPriorFID,'%d %d %d\n',i,j,sign(SMatrix(i,j)));
                end
            end
        end
    end 
    fclose(metsSoftPriorFID);
    metsHardPriorFID=fopen(['metsHardPrior' suffix '.txt'],'w');
    for i=1:size(SMatrix,1)
        for j=1:size(SMatrix,2)
            writeValue=1;
            if(stoichioSoftPrior)
                writeValue=SMatrix(i,j)~=0;
            end
            if(j==size(SMatrix,2))
                fprintf(metsHardPriorFID,'%d\n',writeValue);
            else
                fprintf(metsHardPriorFID,'%d ',writeValue);
            end
        end
    end 
    fclose(metsHardPriorFID);
    
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