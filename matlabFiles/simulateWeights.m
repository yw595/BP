function [T,Y]=simulateWeights(filename)
    WFile=fopen(filename);
    line=fgetl(WFile);
    W=[];
    linenum=1;
    while(line~=-1)
        words=strsplit(line);
        words=words(2:end);
        for i=1:length(words)
            W(linenum,i)=str2num(words{i});
        end
        line=fgetl(WFile);
        linenum=linenum+1;
    end
    W
    [T,Y]=ode15s(@diffEq,[0 1000],5*ones(size(W,1),1));
    
    function dy = diffEq(t,y)
        dy=zeros(size(W,1),1);
        for i=1:length(dy)
            dy(i)=tanh(W(i,:)*y)-y(i);
        end
    end
end