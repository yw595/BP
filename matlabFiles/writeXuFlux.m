[junk1 junk2 totalData]=xlsread('../XuFlux/nchembio.941-S2.xlsx');

names=totalData(3:86,1);
for i=1:length(names)
    names{i}=strrep(names{i},' ','');
end
dataNorm=[];
zeroTimePointIdxs=[2,3,20,21,38,39,56,57];
for i=3:86
    for j=1:length(zeroTimePointIdxs)
        dataNorm(i-2,j)=totalData{i,zeroTimePointIdxs(j)};
    end
end
stdevsNorm=std(dataNorm,0,2);
dataNorm=mean(dataNorm,2);

dataTemp=[];
timeCourseIdxs=[4:19,22:37,40:55,58:73];
for i=3:86
    for j=1:length(timeCourseIdxs)
        dataTemp(i-2,j)=totalData{i,timeCourseIdxs(j)};
    end
end

data=[];stdevs=[];
for i=1:size(dataTemp,1)
    for j=1:size(dataTemp,2)
        if(mod(j,2)==0)
            data(i,j/2)=mean(dataTemp(i,j-1:j));
            stdevs(i,j/2)=std(dataTemp(i,j-1:j),0,2);
        end
    end
end
            

suffix='XuFlux';
minDataVal=min(min(data));
minDataVal=min([dataNorm; minDataVal]);
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