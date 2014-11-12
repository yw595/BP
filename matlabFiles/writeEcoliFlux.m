%cd C:\MinGW\msys\1.0\home\Yiping' Wang'\
initCobraToolbox;
modelEcoli=readCbModel('../EcoliFlux/iJR904.xml');
model=modelEcoli;

%numWrittenRxns is number of rxns to actually write to data and thus infer
numWrittenRxns=length(model.rxns);
%solBiomass is baseline
solBiomass=optimizeCbModel(model,'max','one');
fluxBiomass=solBiomass.x;
fluxDists=[];
optRxnIdxs=1:length(model.rxns);
for i=1:length(optRxnIdxs)
    i
    model=changeObjective(model,{model.rxns{optRxnIdxs(i)}},[1]);
    solIth=optimizeCbModel(model,'max','one');
    fluxDists(:,i)=solIth.x;
end

suffix='EcoliFlux';
epsFlux=.01;

dataFID=fopen(['data' suffix '.txt'],'w');
for i=1:numWrittenRxns%size(fluxDists,1)
    for j=1:size(fluxDists,2)-1
        writeValue=log2((fluxDists(i,j)+epsFlux)/(fluxBiomass(i)+epsFlux));
        if(writeValue<=-9)
            writeValue=-9;
        end
        if(writeValue>=99)
            writeValue=99;
        end
        fprintf(dataFID,'%2.2f ',writeValue);
    end
    if(isempty(j))
        j=1;
    else
        j=j+1;
    end
    writeValue=log2((fluxDists(i,j)+epsFlux)/(fluxBiomass(i)+epsFlux));
    if(writeValue<=-9)
        writeValue=-9;
    end
    if(writeValue>=99)
        writeValue=99;
    end
    fprintf(dataFID,'%2.2f\n',writeValue);
end
fclose(dataFID);

pertFID=fopen(['pert' suffix '.txt'],'w');
for i=1:numWrittenRxns%size(fluxDists,1)
    for j=1:size(fluxDists,2)-1
        fprintf(pertFID,'%2.2f ',0);
    end
    %insurance in case j does not increment on last loop
    if(isempty(j))
        j=1;
    else
        j=j+1;
    end
    if(i==size(fluxDists,1))
        fprintf(pertFID,'%2.2f\n',0);
    else
        fprintf(pertFID,'%2.2f\n',0);
    end
end
fclose(pertFID);

nameFID=fopen(['name' suffix '.txt'],'w');
for i=1:numWrittenRxns%size(fluxDists,1)
    if(i==10)%size(fluxDists,1))
        fprintf(nameFID,'1 %s\n',model.rxns{i});
    else
        fprintf(nameFID,'1 %s\n',model.rxns{i});
    end
end
fclose(nameFID);

priorFID=fopen(['prior' suffix '.txt'],'w');
for i=1:size(fluxDists,1)
end
fclose(priorFID);

inputFID=fopen(['input' suffix '.txt'],'w');
%Session ID
fprintf(inputFID,'%s\n',suffix);
%Nexpts
fprintf(inputFID,'%d\n',length(optRxnIdxs));
%Nnodes
fprintf(inputFID,'%d\n',numWrittenRxns);
%Nwvals
fprintf(inputFID,'%d\n',11);
%MaxWvals
fprintf(inputFID,'%s\n','1.00');
%thresh
fprintf(inputFID,'%s\n','1.00E-06');
%lambda
fprintf(inputFID,'%s\n','3.00');
%beta
fprintf(inputFID,'%s\n','2.00');
%Nprior
fprintf(inputFID,'%d\n',0);
%Nobs
fprintf(inputFID,'%d\n',numWrittenRxns);
%Ndeci
fprintf(inputFID,'%d\n',1);
%opt_maxSkip
fprintf(inputFID,'%d\n',10);
fclose(inputFID);