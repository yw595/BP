function renameInputFiles(suffix,useMets)

system(['rm ' '../' suffix '/name.txt']);
system(['rm ' '../' suffix '/pert.txt']);
system(['rm ' '../' suffix '/prior.txt']);
system(['rm ' '../' suffix '/data.txt']);
%system(['rm ' '../' suffix '/input.txt'])

copyfile(['../' suffix '/fluxName.txt'],['../' suffix '/name.txt']);
copyfile(['../' suffix '/fluxPert.txt'],['../' suffix '/pert.txt']);
if(~exist('useMets','var') && useMets)
    copyfile(['../' suffix '/metsHardPrior.txt'],['../' suffix '/hardPrior.txt']);
    copyfile(['../' suffix '/metsSoftPrior.txt'],['../' suffix '/prior.txt']);
    copyfile(['../' suffix '/metsData.txt'],['../' suffix '/data.txt']);
else
    copyfile(['../' suffix '/fluxHardPrior.txt'],['../' suffix '/hardPrior.txt']);
    copyfile(['../' suffix '/fluxSoftPrior.txt'],['../' suffix '/prior.txt']);
    copyfile(['../' suffix '/fluxData.txt'],['../' suffix '/data.txt']);
end
copyfile(['../' suffix '/input.txt'],'../input.txt');

end