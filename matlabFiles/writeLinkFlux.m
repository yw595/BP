[junk1 junk2 totalData]=xlsread('../LinkFlux/nbt.2489-S2.xlsx');

names=totalData(5:2:117,2);
for i=1:length(names)
    names{i}=strrep(names{i},' ','');
end
dataNorm=[];
for i=5:2:105
    dataNorm(end+1)=totalData{i,17};
end
data=[];
for j=3:9
    for i=5:2:117
        data(i-4,j-2)=totalData{i,j};
    end
end
for j=10:16
    for i=5:2:105
        data(i-4,j-9)=totalData{i,j};
    end
end
stdevsNorm=[];
for i=6:2:106
    stdevsNorm(i-5)=totalData{i,17};
end
stdevs=[];
for j=3:16
    for i=6:2:106
        stdevs(i-5,j-2)=totalData{i,j};
    end
end

names=names(strcmp(totalData(5:2:117,1),'A'));
dataNorm=dataNorm(strcmp(totalData(5:2:105,1),'A'));
data=data(strcmp(totalData(5:2:117,1),'A'));
stdevsNorm=stdevsNorm(strcmp(totalData(5:2:105,1),'A'));
stdevs=stdevs(strcmp(totalData(5:2:105,1),'A'));
data=data(1:end-6,:);
names=names(1:end-6);

suffix='LinkFlux';
minDataVal=min(min(data));
minDataVal=min([dataNorm minDataVal]);
%data=data-minDataVal;
%dataNorm=dataNorm-minDataVal;
epsFlux=.01*min(min(abs(data(data~=0))));

dataFID=fopen(['data' suffix '.txt'],'w');
for i=1:size(data,1)
    for j=1:size(data,2)-1
        writeValue=log2((data(i,j)+epsFlux)/(dataNorm(i)+epsFlux));
        if(writeValue<=-99)
            writeValue=-99;
        end
        if(writeValue>=99)
            writeValue=99;
        end
        fprintf(dataFID,'%2.2f ',writeValue);
    end
    if(isempty(j))
        j=1;
    else
        j=j+1;
    end
    writeValue=log2((data(i,j)+epsFlux)/(dataNorm(i)+epsFlux));
    if(writeValue<=-99)
        writeValue=-99;
    end
    if(writeValue>=99)
        writeValue=99;
    end
    fprintf(dataFID,'%2.2f\n',writeValue);
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
fprintf(inputFID,'%s\n',suffix);
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