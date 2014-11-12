[metsNames fluxNames]=renameMetsAndRxns2();

%parse equations written out in rxnEqs.txt
SMatrix=formSMatrix(metsNames,fluxNames);

rxnInfluenceMatrix=formInfluenceMatrix(fluxNames,SMatrix,metsNames);

%write out storeRxnInfluence as .sif file
realNetworkSif=fopen('realNetwork.sif','w');
for i=1:length(fluxNames)
    for j=1:length(fluxNames)
        if(rxnInfluenceMatrix(i,j)==1)
            fprintf(realNetworkSif,[fluxNames{i} ' activates ' fluxNames{j} '\n']);
        elseif(rxnInfluenceMatrix(i,j)==-1)
            fprintf(realNetworkSif,[fluxNames{i} ' inhibits ' fluxNames{j} '\n']);
        end
    end
end
fclose(realNetworkSif);

%write out storeRxnInfluence as edge attributes file
realAverageEdges=fopen('realAverageEdges.txt','w');
fprintf(realAverageEdges,'EdgeProbability\n');
for i=1:length(fluxNames)
    for j=1:length(fluxNames)
        if(rxnInfluenceMatrix(i,j)==1)
            fprintf(realAverageEdges,[fluxNames{i} ' (activates) ' fluxNames{j} ' = 0.99\n']);
        elseif(rxnInfluenceMatrix(i,j)==-1)
            fprintf(realAverageEdges,[fluxNames{i} ' (inhibits) ' fluxNames{j} ' = -0.99\n']);
        end
    end
end
fclose(realAverageEdges);

[metsNames fluxNames]=renameMetsAndRxns2();

%write out storeRxnInfluence as .sif file
SNetworkSif=fopen('S.sif','w');
for i=1:length(metsNames)
    for j=1:length(fluxNames)
        if(SMatrix(i,j)>0)
            fprintf(SNetworkSif,[fluxNames{j} ' activates ' metsNames{i} '\n']);
        elseif(SMatrix(i,j)<0)
            fprintf(SNetworkSif,[metsNames{i} ' activates ' fluxNames{j} '\n']);
        end
    end
end
fclose(SNetworkSif);

%write out storeRxnInfluence as edge attributes file
SAverageEdges=fopen('SAverageEdges.txt','w');
fprintf(SAverageEdges,'EdgeProbability\n');
for i=1:length(metsNames)
    for j=1:length(fluxNames)
        if(SMatrix(i,j)>0)
            fprintf(SAverageEdges,[fluxNames{j} ' (activates) ' metsNames{i} ' = 0.99\n']);
        elseif(SMatrix(i,j)<0)
            fprintf(SAverageEdges,[metsNames{i} ' (activates) ' fluxNames{j} ' = -0.99\n']);
        end
    end
end
fclose(SAverageEdges);

%write out storeRxnInfluence as edge attributes file
SAverageNodes=fopen('SAverageNodes.txt','w');
fprintf(SAverageNodes,'Type\n');
for i=1:length(metsNames)
    fprintf(SAverageNodes,[metsNames{i} ' = met\n']);
end
for i=1:length(fluxNames)
    fprintf(SAverageNodes,[fluxNames{i} ' = flux\n']);
end
fclose(SAverageNodes);