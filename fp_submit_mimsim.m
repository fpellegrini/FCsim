

fp_addpath
cd ~/matlab/fp/meg_lfp/
nit = 100;
varyParam = [1:8];%defaults, snr, interactions, everything else 

for ip =7 %1,2,3,4,5,6, 7 
    mgsub({},@fp_eval_mim_struct_sim,{ip},'qsub_opts',['-l h_vmem=16G -t 1-' num2str(nit) ]) %' -q 2jobs'
    pause(60)
end