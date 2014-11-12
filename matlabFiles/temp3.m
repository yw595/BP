%inputFile1='StoichioPrior2Again\Model_1.sif';
inputFile=fopen('StoichioPrior2again\Model_1edges.txt');

edges={};
influenceStrengths=[];
fgetl(inputFile);
line=fgetl(inputFile);
while(line~=-1)
    firstNode=line(1:regexp(line,'(')-2);
    secondNode=line(regexp(line,')')+2:regexp(line,'=')-2);
    influence=line(regexp(line,'(')+1:regexp(line,')')-1);
    influenceStrength=str2num(line(regexp(line,'=')+3:end));
    edges{end+1,1}=firstNode;
    edges{end,2}=secondNode;
    edges{end,3}=influence;
    influenceStrengths(end+1)=influenceStrength;
    line=fgetl(inputFile);
end
fclose(inputFile);

[junk sortIdxs]=sort(abs(influenceStrengths));
influenceStrengths=influenceStrengths(sortIdxs(end-99:end));
edges=edges(sortIdxs(end-99:end),:);

outputFile=fopen('StoichioPrior2again\Modelfiltered_1.sif','w');
for i=1:length(edges)
    fprintf(outputFile,'%s %s %s\n',edges{i,1},edges{i,3},edges{i,2});
end
fclose(outputFile);

outputFile=fopen('StoichioPrior2again\Modelfiltered_1edges.txt','w');
fprintf(outputFile,'EdgeProbability\n');
for i=1:length(edges)
    fprintf(outputFile,'%s (%s) %s =  %f\n',edges{i,1},edges{i,3},edges{i,2},influenceStrengths(i));
end
fclose(outputFile);