function rxnInfluenceMatrix=formRxnInfluenceMatrix(metsNames,rxnNames,SMatrix)
    
    %form rxnInfluenceMatrix by dot-product of S-matrix columns for reaction,
    %filter out elements at CO2 rows, and remove Fba and Tpi interactions since
    %are ambiguous
    rxnInfluenceMatrix=zeros(length(rxnNames),length(rxnNames));
    for i=1:length(rxnNames)
        for j=1:length(rxnNames)
            if(i~=j)
                interactionVector=SMatrix(:,i).*SMatrix(:,j);
                if(size(SMatrix,1)>=find(strcmp(metsNames,'CO2')))
                    interactionVector=interactionVector-...
                    SMatrix(strcmp(metsNames,'CO2'),i).*SMatrix(strcmp(metsNames,'CO2'),j)-...
                    SMatrix(strcmp(metsNames,'ATP'),i).*SMatrix(strcmp(metsNames,'ATP'),j)-...
                    SMatrix(strcmp(metsNames,'ADP'),i).*SMatrix(strcmp(metsNames,'ADP'),j);
                end
                if((sum(strcmp(rxnNames{i},'Tpi'))==0 && sum(strcmp(rxnNames{j},'Fba'))==0) ...
                && (sum(strcmp(rxnNames{i},'Fba'))==0 && sum(strcmp(rxnNames{j},'Tpi'))==0) ...
                && (sum(strcmp(rxnNames{i},'PykA'))==0 && sum(strcmp(rxnNames{j},'Pgk;Pgm;Eno'))==0) ...
                && (sum(strcmp(rxnNames{i},'Pgk;Pgm;Eno'))==0 && sum(strcmp(rxnNames{j},'PykA'))==0))
                    %check if metabolite interactions are of different sign,
                    %should only be Tpi and Fba
                    if(sum(interactionVector<0)>0 && sum(interactionVector>0)>0)
                        disp('AMBIGUOUS');
                        disp(rxnNames{i});
                        disp(rxnNames{j});
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