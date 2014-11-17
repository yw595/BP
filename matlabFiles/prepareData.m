function data = prepareData(data,stds,duplicate,localMin,globalMin,useLog)  
    %use standard deviations to replicate data 5 times with Gaussian noise
    dataTemp=[];
    for i=1:size(data,1)
        for j=1:size(data,2)
            for k=1:5
                dataTemp(i,(j-1)*5+k)=data(i,j)+stds(i,j)*randn(1,1);
            end
        end
    end
    if(duplicate)
        data=dataTemp;
    end

    %calculate minimum absolute data value, subtract from all data so all
    %data are positive, using either localMin for each row individually, or
    %globalMin for all data at once
    if(localMin)
        for i=1:size(data,1)
            minDataVal=min(data(i,:));
            data(i,:)=data(i,:)-minDataVal;
        end
    end
    if(globalMin)
        minDataVal=min(min(data));
        data=data-minDataVal;
    end
    
    %log transform data if useLog, taking first column as baseline and
    %discarding it
    epsVal=min(min(abs(data(data~=0))));
    if(useLog)
        dataNorm=data(:,1);
        if(duplicate)
            data=data(:,6:end);
        else
            data=data(:,2:end);
        end
        for i=1:size(data,2)
            for j=1:size(data,1)
                if(data(j,i)==0 || dataNorm(j)==0)
                    data(j,i)=log((data(j,i)+epsVal)/(dataNorm(j)+epsVal));
                else
                    data(j,i)=log(data(j,i)/dataNorm(j));
                end
            end
        end
    end
end