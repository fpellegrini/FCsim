function fp_submit_mimsim

fp_addpath
cd ~/matlab/fp/
nit = 100;
varyParam = [1:7];%defaults, snr, interactions, everything else 

for ip =1
    mgsub({},@fp_eval_mim_struct_sim,{ip},'qsub_opts',['-l h_vmem=16G -t 1-' num2str(nit) ]) %' -q 2jobs'
    pause(2)
end