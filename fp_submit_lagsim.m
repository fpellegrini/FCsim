

fp_addpath
cd ~/matlab/fp/meg_lfp/
nit = 100;
varyParam = [1 7];%defaults, snr, interactions, everything else 

for ip = 1 %todo: ip7
    mgsub({},@fp_eval_lagsim,{ip},'qsub_opts',['-l h_vmem=16G -t 1-' num2str(nit) ]) %' -q 2jobs'
    pause(60)
end