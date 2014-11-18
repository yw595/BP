function matrix = readMatrix(inputFile,delimiter)
    matrix = [];
    inputFID=fopen([suffix '\' inputFile]);
    line=fgetl(inputFID);
    i=0;
    i=i+1;
    while(line~=-1)
        words=strsplit(line,delimiter);
        words=words(2:end);
        for j=1:length(words)
            matrix(i,j)=str2num(words{j});
        end
        line=fgetl(modelFID);
        i=i+1;
    end
    fclose(inputFID);
end