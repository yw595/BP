function rxnInfluenceMatrix=formInfluenceMatrix(fluxNames,SMatrix,metsNames)
    %form rxnInfluenceMatrix by dot-product of S-matrix columns for reaction,
    %filter out elements at CO2 rows, and remove Fba and Tpi interactions since
    %are ambiguous
    rxnInfluenceMatrix=zeros(length(fluxNames),length(fluxNames));
    for i=1:length(fluxNames)
        for j=1:length(fluxNames)
            if(i~=j)
                interactionVector=SMatrix(:,i).*SMatrix(:,j);
                if(size(SMatrix,1)>=find(strcmp(metsNames,'CO2')))
                    interactionVector=interactionVector-...
                    SMatrix(strcmp(metsNames,'CO2'),i).*SMatrix(strcmp(metsNames,'CO2'),j)-...
                    SMatrix(strcmp(metsNames,'ATP'),i).*SMatrix(strcmp(metsNames,'ATP'),j)-...
                    SMatrix(strcmp(metsNames,'ADP'),i).*SMatrix(strcmp(metsNames,'ADP'),j);
                end
                if((sum(strcmp(fluxNames{i},'Tpi'))==0 && sum(strcmp(fluxNames{j},'Fba'))==0) ...
                && (sum(strcmp(fluxNames{i},'Fba'))==0 && sum(strcmp(fluxNames{j},'Tpi'))==0) ...
                && (sum(strcmp(fluxNames{i},'PykA'))==0 && sum(strcmp(fluxNames{j},'Pgk;Pgm;Eno'))==0) ...
                && (sum(strcmp(fluxNames{i},'Pgk;Pgm;Eno'))==0 && sum(strcmp(fluxNames{j},'PykA'))==0))
                    %check if metabolite interactions are of different sign,
                    %should only be Tpi and Fba
                    if(sum(interactionVector<0)>0 && sum(interactionVector>0)>0)
                        disp('AMBIGUOUS');
                        disp(fluxNames{i});
                        disp(fluxNames{j});
                    elseif(sum(interactionVector)==1)
                        rxnInfluenceMatrix(i,j)=-1;
                    elseif(sum(interactionVector)==-1)
                        rxnInfluenceMatrix(i,j)=1;
                    end
                end
            end
        end
    end
end