function extendedInfluenceMatrices = makeExtendedInfluenceMatrices(originalInfluenceMatrix,maxZ)

    extendedInfluenceMatrices = {};
    absInfluenceMatrix = abs(originalInfluenceMatrix);
    extendedInfluenceMatrices{1} = originalInfluenceMatrix;
    for z = 2:maxZ
        influenceMatrix = extendedInfluenceMatrices{z};
        for i=1:size(originalInfluenceMatrix,1)
            for j=1:size(originalInfluenceMatrix,2)
                [junk1,path,junk2]=graphshortestpath(sparse(absInfluenceMatrix),i,j);
                if(length(path)<=cutoff+1 && length(path)>1)
                    sign=1;
                    for k=1:length(path)-1
                        sign=sign*originalInfluenceMatrix(path(k),path(k+1));
                    end
                    influenceMatrix(i,j)=sign;
                end
            end
        end
        extendedInfluenceMatrices{z} = influenceMatrix;
    end
end