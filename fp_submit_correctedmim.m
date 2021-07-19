
fp_addpath
cd ~/matlab/fp/meg_lfp/
nit = 100;
varyParam = 1

mgsub({},@fp_add_corrected_pips,{},'qsub_opts',['-l h_vmem=16G -t 1-' num2str(nit) ]) %' -q 2jobs'
