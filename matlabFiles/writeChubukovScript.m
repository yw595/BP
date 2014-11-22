useMets=0;
runNames={'ChubukovTest'};
optionsArrays={[1 0 0 1 0 0 0 1 1 1 2 0 0 1 10]};

if(0)
runNames={};
optionsArrays={};
for i=1:10
    for j=1:10
        runNames{end+1}=['testCostFluxi' num2str(i) 'j' num2str(j)];
        optionsArrays{end+1}=[0 1 0 1 0 0 0 0 1 1 2 0 0 i j];
    end
end
end

stats={};
for i=1:length(runNames)
    optionsArray=optionsArrays{i};
    localMin=optionsArray(1);globalMin=optionsArray(2);perturbations=optionsArray(3);duplicate=optionsArray(4);
    influenceHardPrior=optionsArray(5);correlationHardPrior=optionsArray(6);influenceSoftPrior=optionsArray(7);
    correlationSoftPrior=optionsArray(8);useFluxLog=optionsArray(9);useMetsLog=optionsArray(10);
    deciOption=optionsArray(11);stoichioSoftPrior=optionsArray(12);useLinearAnsatz=optionsArray(13);
    beta=optionsArray(14);lambda=optionsArray(15);
    
    suffix=runNames{i};
    writeChubukov(suffix,[],[32 35],localMin,globalMin,perturbations,duplicate,influenceHardPrior,...
        correlationHardPrior,influenceSoftPrior,correlationSoftPrior,useFluxLog,useMetsLog,deciOption,...
        stoichioSoftPrior,useLinearAnsatz,beta,lambda);
    
    renameInputFiles(suffix,useMets);
    eval('cd ..');
    system('make');
    disp(['.' sprintf('\\doBP_FULL\n0\n%s\n', suffix)])
    system(['.' sprintf('\\doBP_FULL\n')]);
    eval('cd matlabFiles');
    break;
    
    if(0)
    [signedRecall signedPrecision signedCounts unsignedRecall unsignedPrecision unsignedCounts] ...
        =validateCorrelate([suffix '\Model_1edges.txt'],[],[32 35],useMets);
    stats{i,1}=signedRecall;
    stats{i,2}=signedPrecision;
    stats{i,3}=signedCounts;
    stats{i,4}=unsignedRecall;
    stats{i,5}=unsignedPrecision;
    stats{i,6}=unsignedCounts;
    stats{i,7}=(i-mod(i,10))/10;
    stats{i,8}=mod(i,10);
    if(stats{i,8}==0)
        stats{i,8}=10;
    end
    stats(i,9)=i;
    end
end