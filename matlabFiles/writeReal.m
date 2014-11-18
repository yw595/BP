[metsNames fluxNames]=renameMetsAndRxns('includeExtraMets');

SMatrix=formSMatrix(metsNames,fluxNames);
rxnInfluenceMatrix=formInfluenceMatrix(fluxNames,SMatrix,metsNames);

writeMatrix(rxnInfluenceMatrix,'realNetwork.sif',0,1,0,.25,fluxNames,fluxNames);
realAverageEdges = .99*rxnInfluenceMatrix>0 -.99*rxnInfluenceMatrix<0;
writeMatrix(realAverageEdges,'realAverageEdges.txt',0,0,1,.25,fluxNames,fluxNames);

writeMatrix(SMatrix,'S.sif',0,1,0,.25,metsNames,fluxNames);
SAverageEdges = .99*SMatrix>0 -.99*SMatrix<0;
writeMatrix(SAverageEdges,'SAverageEdges.txt',0,0,1,.25,metsNames,fluxNames);

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