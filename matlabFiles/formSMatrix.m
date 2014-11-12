function SMatrix=formSMatrix(metsNames,rxnNames)
    %parse equations written out in rxnEqs.txt
    SMatrix=zeros(length(metsNames),length(rxnNames));
    rxnEqsFID=fopen('../Chubukov/rxnEqs.txt');
    line=fgetl(rxnEqsFID);
    while(line~=-1)
        rxnName=line(1:regexp(line,':')-1);
        %get reactants and products parts of line
        reactantsLine=line(regexp(line,':')+1:regexp(line,'-')-1);
        productsLine=line(regexp(line,'>'):end);
        reactants=[];
        products=[];
        %if reactant and product parts are not empty, split them on +
        if(length(reactantsLine)~=1)
            reactantsLine=reactantsLine(2:end-1);
            reactants=strsplit(reactantsLine,' + ');
        end
        if(length(productsLine)~=1)
            productsLine=productsLine(3:end);
            products=strsplit(productsLine,' + ');
        end
        %if the reaction Name matches the storeNames with spaces, put one in
        %corresponding position in S matrix
        for i=1:length(products)
            if(sum(strcmp(rxnNames,rxnName))~=0 && sum(strcmp(metsNames,products{i}))~=0)
                SMatrix(strcmp(metsNames,products{i}),strcmp(rxnNames,rxnName))=1;
            end
        end
        for i=1:length(reactants)
            if(sum(strcmp(rxnNames,rxnName))~=0 && sum(strcmp(metsNames,reactants{i}))~=0)
                SMatrix(strcmp(metsNames,reactants{i}),strcmp(rxnNames,rxnName))=-1;
            end
        end
        line=fgetl(rxnEqsFID);
    end
    fclose(rxnEqsFID);
end