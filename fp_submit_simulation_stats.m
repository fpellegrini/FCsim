fp_addpath
cd ~/matlab/fp/meg_lfp/
nit = 100;
mgsub({},@fp_simulation_stats2,{},'qsub_opts',['-l h_vmem=16G -t 1-' num2str(nit) ]) %' -q 2jobs'
