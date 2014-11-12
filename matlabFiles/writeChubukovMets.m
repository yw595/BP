suffix='ChubukovMets';

[junk1 junk2 totalData]=xlsread('../ChubukovFlux/inline-supplementary-material-2.xlsx',3);

CBNames={'arginine','asparagine','glutamine','disaccharides','tyrosine','adenine',...
    'phenylalanine','aspartate','glutamate','glucose-6-phosphate (G6P)',...
    'ribose-5-phosphate (R5P)','glycerol-3-phosphate (G3P)','ribulose-5-phosphate (Ru5P)',...
    'xylulose-5-phosphate (Xu5P)','dihydroxyacetone-phosphate (DHAP)','ribose-1-phosphate (R1P)',...
    'guanosine monophosphate (GMP)','pantothenate','adenosine monophosphate (AMP)','succinate',...
    'malate','alpha-ketoglutarate','fumarate','phospho-glycerate (2PG+3PG)','6-phospho-gluconate',...
    'guanosine diphosphate (GDP)','adenosine diphosphate (ADP)','phosphoenolpyruvate (PEP)','fructose-bis-phosphate (FBP)',...
    'citrate','isocitrate','guanosine triphosphate (GTP)','adenosine triphosphate (ATP)','bis-phosho-glycerate (BPG)'};
NewNames={'','','','','','',...
    '','','GLUT','G6P',...
    'R5P','','Ru5P',...
    'X5P','DHAP','',...
    '','','','',...
    'MAL','AKG','','','6PGN',...
    '','','PEP','FBP',...
    'CITR','','','','13PG'};
CBNamesToNewNames=containers.Map(CBNames,NewNames);

names=totalData(5:38,1);
for i=1:length(names)
    if(~strcmp(CBNamesToNewNames(names{i}),''))
        names{i}=CBNamesToNewNames(names{i});
    end
    names{i}=strrep(names{i},' ','');
end

data=[];
for j=2:9
    for i=5:38
        if(strcmp(totalData{i,j},'n/a'))
            data(i-4,j-1)=0;
        else
            data(i-4,j-1)=totalData{i,j};
        end
    end
end

dataFID=fopen(['mets' suffix '.txt'],'w');
for i=1:size(data,1)
    for j=1:size(data,2)
        if(j==size(data,2))
            fprintf(dataFID,'%2.2f\n',data(i,j));
        else
            fprintf(dataFID,'%2.2f ',data(i,j));
        end
    end
end
fclose(dataFID);

%parse equations written out in rxnEqs.txt
dataFID=fopen(['metLinks' suffix '.txt'],'w');
SMatrix2=[];
for i=1:size(SMatrix,1)
    metIdx=strcmp(NewNames,metabolites{i});
    SMatrix2(metIdx,:)=SMatrix(i,:);
end
for i=1:size(SMatrix2,1)
    for j=1:size(SMatrix2,2)
        if(j==size(SMatrix2,2))
            fprintf(dataFID,'%d\n',SMatrix2(i,j));
        else
            fprintf(dataFID,'%d ',SMatrix2(i,j));
        end
    end
end

line=fgetl(rxnEqsFID);
while(line~=-1)
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
    writeArray=zeros(1,13);
    for i=1:length(products)
        writeArray(strcmp(NewNamesShort,products{i}))=1;
    end
    for i=1:length(reactants)
        writeArray(strcmp(NewNamesShort,reactants{i}))=-1;
    end
    for i=1:length(writeArray)-1
        fprintf(dataFID,'%d ',writeArray(i));
    end
    fprintf(dataFID,'%d\n',writeArray(i));
    line=fgetl(rxnEqsFID);
end
fclose(dataFID);