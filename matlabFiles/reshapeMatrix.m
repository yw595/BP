function matrix = reshapeMatrix(oldMatrix)
    matrix = zeros(size(oldMatrix,1)*size(oldMatrix,2),3);
    for i=1:size(oldMatrix,1)
        for j=1:size(oldMatrix,2)
            matrix((i-1)*size(oldMatrix,2)+j,:)=[i j oldMatrix(i,j)];
        end
    end
    matrix( matrix(:,3)==0,:)=[];
end