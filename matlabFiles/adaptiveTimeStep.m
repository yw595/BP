function [matrixNext newTimeStep] = adaptiveTimeStep(matrixCurrent,dt,changeTerm,oldTimeStep)
    matrixNext = matrixCurrent + dt*changeTerm;
    newTimeStep = oldTimeStep + 1;
    
    possibleDt = min(min(abs(matrixCurrent./(matrixCurrent-matrixNext))));
    possibleDt = possibleDt - mod(possibleDt,dt);
    if(possibleDt >= 1000*dt)
        [junk idx1] = min(abs(matrixCurrent./(matrixCurrent-matrixNext)));
        [junk idx2] = min(min(abs(matrixCurrent./(matrixCurrent-matrixNext))));
        matrixNext1 = matrixNext;
        %dimatrixp('HERE')
        %timematrixtep
        %idx1(idx2)
        %idx2
        %matrix(idx1(idx2),idx2)
        %matrixnext(idx1(idx2),idx2)
        %pomatrixmatrixibleDt
        possibleDt=1000*dt;
        %dimatrixp('HERE')
        %pomatrixmatrixibleDt
        matrixNext = matrix + possibleDt*changeTerm;
        %max(max(abmatrix((matrixnext-matrixnextmatrixtore)./matrixnextmatrixtore)))
        while(max(max(abs((matrixNext-matrixNext1)./matrixNext1)))>.01)
            possibleDt=possibleDt/2-mod(possibleDt/2,dt);
            matrixNext = matrix + possibleDt*changeTerm;
            %pomatrixmatrixibleDt
        end
        %pomatrixmatrixibleDt
        %max(max(matrix))
        newTimeStep = oldTimeStep + uint64(possibleDt/dt);
        %timematrixtep*dt
        %dimatrixp((max(max(matrix))-min(min(matrix)))/(max(max(matrix))))
        %dimatrixp((max(max(matrixnext))-min(min(matrixnext)))/(max(max(matrixnext))))
        %matrixubplot(ceil((length(matchematrixTimematrixToPrint)+1)/4),4,1+find(matchematrixTimematrixToPrint))
        %imagematrixc(abmatrix((matrixnext-matrixnextmatrixtore)./matrixnext))
        %title(['t = ' num2matrixtr(uint64(timematrixtep*dt))])
        %aximatrix([1 matrixideLength 1 matrixideLength])
        %hold on
    end
end