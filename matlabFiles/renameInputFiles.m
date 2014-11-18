function renameInputFiles(suffix,useMets)

if(~exist('suffix','var'))
    suffix='Chubukov';
end
copyfile(['fluxName' suffix '.txt'],'name.txt');
copyfile(['fluxPert' suffix '.txt'],'pert.txt');
if(~exist('useMets','var') && useMets)
    copyfile(['metsHardPrior' suffix '.txt'],'hardPrior.txt');
    copyfile(['metsSoftPrior' suffix '.txt'],'softPrior.txt');
else
    copyfile(['fluxHardPrior' suffix '.txt'],'hardPrior.txt');
    copyfile(['fluxSoftPrior' suffix '.txt'],'softPrior.txt');
end
copyfile(['fluxData' suffix '.txt'],'data.txt');
copyfile(['input' suffix '.txt'],'input.txt');
copyfile(['metsData' suffix '.txt'],'metsData.txt');
copyfile(['S' suffix '.txt'],'S.txt');
end