%submit to nsg cluster

addpath(genpath('./'))
nit = 100;
seeds = (1:nit)*17; 

for iit = 1:nit
    fprintf(['starting job ' num2str(iit) '\n'])
    job{iit} = batch(@fp_simulation_stats2,0,{seeds(iit),iit},'Pool',1);
end

for iit =1:nit
    fprintf(['waiting for job ' num2str(iit) '\n'])
    wait(job{nit})
end
