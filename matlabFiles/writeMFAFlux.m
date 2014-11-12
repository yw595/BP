suffixes={'1OnlyRedZeroBGMin1PercentHigh','1OnlyRedZeroBGMin1PercentLow'};
names={};
data=[];
dataNorm=[];
for z=1:2
    suffix=suffixes{z};
    modelReactionsToOutputReactions=containers.Map;
    modelReactions={};
    modelReactionsReversibility={};
    
    %read in model reactions, reversibility for each one, and map to output
    %reactions
    modelFID=fopen(['../MFAFlux/model' suffix '.txt']);
    line=fgetl(modelFID);
    while(line~=-1)
        if(sum(regexp(line,'^R'))~=0)
            words=strsplit(line,'\t');
            
            reaction=words{1};
            modelReactions{end+1}=reaction;
            modelReactionsReversibility{end+1}=words{4};
            outputReaction=reaction(2:end);
            if(regexp(outputReaction,'FR$'))
                outputReaction=outputReaction(1:end-2);
            elseif(regexp(outputReaction,'R$') | regexp(outputReaction,'F$'))
                outputReaction=outputReaction(1:end-1);
            end
            modelReactionsToOutputReactions(reaction)=outputReaction;
        end
        line=fgetl(modelFID);
    end
    fclose(modelFID);
    
    %read in corresponding results_PE file, map net fluxes, subtracting reverse, for each diagramreaction
    modelFID=fopen(['../MFAFlux/results_PE' suffix '.txt']);
    line=fgetl(modelFID);
    simFluxes=[];
    fluxes=0;
    prevwords={};
    dataAndNamesIdx=1;
    while(line~=-1)
        
        if(sum(regexp(line,'estimated fluxes')))
            fluxes=1;
            line=fgetl(modelFID);
            continue;
        end
        if(sum(regexp(line,'chi-square cut off')))
            fluxes=0;
        end
        
        %parse reactionNum for first word, use to index into
        %modelReactionsReversibility. If forward, we may simply append flux to
        %all diagramreactions and write to file. If R, use prevwords from previous line to
        %and subtract reverse fluxes, then append and write flux.
        if(fluxes)
            words=strsplit(line,'\t');
            reactionNum=words{1};
            reversibility=modelReactionsReversibility{str2num(reactionNum(2:end))};
            if(strcmp(reversibility,'F'))
                outputReaction=modelReactionsToOutputReactions(modelReactions{str2num(reactionNum(2:end))});
                if(z==1)
                    dataNorm(dataAndNamesIdx)=str2num(words{2});
                else
                    data(dataAndNamesIdx,z-1)=str2num(words{2});
                end
                names{dataAndNamesIdx}=outputReaction;
                dataAndNamesIdx=dataAndNamesIdx+1;
            elseif(strcmp(reversibility,'R'))
                outputReaction=modelReactionsToOutputReactions(modelReactions{str2num(reactionNum(2:end))});
                if(z==1)
                    dataNorm(dataAndNamesIdx)=str2num(prevwords{2})-str2num(words{2});
                else
                    data(dataAndNamesIdx,z-1)=str2num(prevwords{2})-str2num(words{2});
                end
                names{dataAndNamesIdx}=outputReaction;
                dataAndNamesIdx=dataAndNamesIdx+1;
            end
            prevwords=words;
        end
        line=fgetl(modelFID);
    end
    fclose(modelFID);
end

suffix='MFAFlux';
epsFlux=.01*min(min(abs(data(data~=0))));

dataFID=fopen(['data' suffix '.txt'],'w');
for i=1:size(data,1)
    for j=1:size(data,2)-1
        fprintf(dataFID,'%2.2f ',log2((data(i,j)+epsFlux)/(dataNorm(i)+epsFlux)));
    end
    if(isempty(j))
        j=1;
    else
        j=j+1;
    end
    if(i==size(data,1))
        fprintf(dataFID,'%2.2f\n',log2((data(i,j)+epsFlux)/(dataNorm(i)+epsFlux)));
    else
        fprintf(dataFID,'%2.2f\n',log2((data(i,j)+epsFlux)/(dataNorm(i)+epsFlux)));
    end
end
fclose(dataFID);

pertFID=fopen(['pert' suffix '.txt'],'w');
for i=1:size(data,1)
    for j=1:size(data,2)-1
        fprintf(pertFID,'%2.2f ',0);
    end
    if(isempty(j))
        j=1;
    else
        j=j+1;
    end
    if(i==size(data,1))
        fprintf(pertFID,'%2.2f\n',0);
    else
        fprintf(pertFID,'%2.2f\n',0);
    end
end
fclose(pertFID);

nameFID=fopen(['name' suffix '.txt'],'w');
for i=1:size(data,1)
    if(i==size(data,1))
        fprintf(nameFID,'1 %s\n',names{i});
    else
        fprintf(nameFID,'1 %s\n',names{i});
    end
end
fclose(nameFID);

priorFID=fopen(['prior' suffix '.txt'],'w');
for i=1:size(data,1)
end
fclose(priorFID);

inputFID=fopen(['input' suffix '.txt'],'w');
fprintf(inputFID,'%s\n','FirstTry');
fprintf(inputFID,'%d\n',size(data,2));
fprintf(inputFID,'%d\n',size(data,1));
fprintf(inputFID,'%d\n',11);
fprintf(inputFID,'%s\n','1.00');
fprintf(inputFID,'%s\n','1.00E-06');
fprintf(inputFID,'%s\n','3.00');
fprintf(inputFID,'%s\n','2.00');
fprintf(inputFID,'%d\n',0);
fprintf(inputFID,'%d\n',size(data,1));
fprintf(inputFID,'%d\n',1);
fprintf(inputFID,'%d\n',10);
fclose(inputFID);