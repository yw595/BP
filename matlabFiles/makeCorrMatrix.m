function corrMatrix = makeCorrMatrix(dataMatrix1,dataMatrix2)

    if(~exist('dataMatrix2','var'))
        dataMatrix2 = dataMatrix1;
    end
    corrMatrix=zeros(size(dataMatrix1,1),size(dataMatrix2,1));
    for i=1:size(dataMatrix1,1)
        for j=1:size(dataMatrix2,1)
            row1=dataMatrix1(i,dataMatrix1(i,:)~=0 & dataMatrix2(j,:)~=0 );
            row2=dataMatrix2(j,dataMatrix1(i,:)~=0 & dataMatrix2(j,:)~=0 );
            if(length(row1)>3)
                tempCorr=corr([row1' row2'],'type','Spearman');
                if(~all(row1==row1(1)) && ~all(row2==row2(1)))
                    corrMatrix(i,j)=tempCorr(1,2);
                end
            end
        end
    end
end