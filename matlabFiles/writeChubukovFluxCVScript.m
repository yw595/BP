system('rm -r ChubukovFluxCV');
system('mkdir ChubukovFluxCV');
for i=1:7
    suffix=['ChubukovFluxCV' num2str(i)];
    writeChubukovFlux(suffix,[i]);
    copyfile(['name' suffix '.txt'],'name.txt');
    copyfile(['pert' suffix '.txt'],'pert.txt');
    copyfile(['prior' suffix '.txt'],'prior.txt');
    copyfile(['data' suffix '.txt'],'data.txt');
    copyfile(['input' suffix '.txt'],'input.txt');
    system('bin\\Debug\\BeliefPropagationCB.exe < input.txt');
    currentDate=clock;
    system(['mv ' num2str(currentDate(1)) '0' num2str(currentDate(2)) ...
        num2str(currentDate(3)) '_' suffix ' ChubukovFluxCV']);
end