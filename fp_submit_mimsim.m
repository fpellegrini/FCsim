

fp_addpath
cd ~/matlab/fp/meg_lfp/
nit = 100; %100 iterations per experiment
varyParam = [1:5 7 9]; %different experimental setups 

for ip = 7
    mgsub({},@fp_eval_mim_struct_sim,{ip},'qsub_opts',['-l h_vmem=16G -t 1-' num2str(nit) ]) %' -q 2jobs'
    pause(60*2)
end