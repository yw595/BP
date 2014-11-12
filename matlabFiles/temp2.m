stats=statsMetsNoPrior;
highestIdx=0;
highestPrecision=0;
for i=1:size(stats,1)
    precisionArray=stats{i,2};
    if(precisionArray(1)>highestPrecision)
        highestIdx=i;
        highestPrecision=precisionArray(1);
    end
end