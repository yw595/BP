function [metsNames rxnNames]=renameMetsAndRxns(includeExtraMets)
    %where possible, use same abbreviations as Chubukov data, else make up
    %own avvrebiations
    %includeExtraMets controls whether unmeasured metabolites in network
    %are included
    metsNames={'','','','','','',...
    '','','GLUT','G6P',...
    'R5P','','Ru5P',...
    'X5P','DHAP','',...
    '','','','SUC',...
    'MAL','AKG','FUM','PG','6PGN',...
    'GDP','ADP','PEP','FBP',...
    'CITR','ISO','GTP','ATP','13PG'};
    metsNames{1}='ARG';metsNames{2}='ARGN';metsNames{3}='GLN';metsNames{4}='DIS';
    metsNames{5}='TYR';metsNames{6}='ADN';metsNames{7}='PHN';metsNames{8}='ASP';
    metsNames{12}='G3P';metsNames{16}='R1P';metsNames{17}='GMP';metsNames{18}='PTN';
    metsNames{19}='AMP';
    if(exist('includeExtraMets','var') && strcmp(includeExtraMets,'includeExtraMets'))
        metsNames{end+1}='F6P';
        metsNames{end+1}='GAP';
        metsNames{end+1}='PYR';
        metsNames{end+1}='E4P';
        metsNames{end+1}='S7P';
        metsNames{end+1}='AcCoA';
        metsNames{end+1}='CO2';
        metsNames{end+1}='OAA';
        metsNames{end+1}='ACE';
    end
    
    rxnNames={};
    rxnEqsFID=fopen('../ChubukovData/rxnEqs.txt');
    line=fgetl(rxnEqsFID);
    while(line~=-1)
        rxnNames{end+1}=line(1:regexp(line,':')-1);
        line=fgetl(rxnEqsFID);
    end
    fclose(rxnEqsFID);
end