function writeMatrix(matrix,outputFile,writeInteger,writeSif,writeEdges,sigThresh,names1,names2)
    outputFID = fopen(outputFile,'w');
    if(exist('writeSif','var') && writeSif)
        for i=1:size(matrix,1)
            for j=1:size(matrix,2)
                if(abs(matrix(i,j))>=sigThresh)
                    if(matrix(i,j)<0)
                        fprintf(outputFID,'%s (inhibits) %s\n',names1{j},names2{i});
                        fprintf(ModelAverageEdgesFID,'%s (inhibits) %s = %f\n',metsNames{j},fluxNames{i},WAverage(i,j));
                    end
                    if(matrix(i,j)>0)
                        fprintf(outputFID,'%s (activates) %s\n',names1{j},names2{i});
                        fprintf(ModelAverageEdgesFID,'%s (activates) %s = %f\n',metsNames{j},fluxNames{i},WAverage(i,j));
                    end
                end
            end
        end
    elseif(exist('writeEdges','var') && writeEdges)
        for i=1:size(matrix,1)
            for j=1:size(matrix,2)
                if(abs(matrix(i,j))>=sigThresh)
                    if(matrix(i,j)<0)
                        fprintf(outputFID,'%s (inhibits) %s = %f\n',names1{j},names2{i},matrix(i,j));
                    end
                    if(matrix(i,j)>0)
                        fprintf(outputFID,'%s (activates) %s = %f\n',names1{j},names2{i},matrix(i,j));
                    end
                end
            end
        end
    else
        for i=1:size(matrix,1)
            for j=1:size(matrix,2)
                if(exist('writeInteger','var') && writeInteger==1)
                    if(j==size(matrix,2))
                        fprintf(outputFID,'%d\n',matrix(i,j));
                    else
                        fprintf(outputFID,'%d ',matrix(i,j));
                    end
                end
                if(j==size(matrix,2))
                    fprintf(outputFID,'%2.2f\n',matrix(i,j));
                else
                    fprintf(outputFID,'%2.2f ',matrix(i,j));
                end
            end
        end
    end
    fclose(outputFID);
end