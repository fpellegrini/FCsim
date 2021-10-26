fp_addpath
cd ~/matlab/fp/meg_lfp/

mgsub({},@fp_addpip,{},'qsub_opts','-l h_vmem=16G')
pause(60)