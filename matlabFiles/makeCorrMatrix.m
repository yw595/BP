function corrMatrix = makeCorrMatrix(dataMatrix1,dataMatrix2)

    if(~exist('dataMatrix2','var'))
        dataMatrix2 = dataMatrix2;
    end
    corrMatrix=zeros(size(dataMatrix1,1),size(dataMatrix2,1));
    for i=1:size(dataMatrix1,1)
        for j=1:size(dataMatrix2,1)
            row1=dataMatrix1(i, dataMatrix1(i,:) & dataMatrix2(j,:) );
            row2=dataMatrix2(j, dataMatrix1(i,:) & dataMatrix2(j,:) );
            if(length(row1)>3)
                tempCorr=corr([row1' row2'],'type','Spearman');
                corrMatrix(i,j)=tempCorr(1,2);
            end
        end
    end
end