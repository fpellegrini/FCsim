function qsubthis(jobname,input, njobs,nsecs)

if nargin<2
    dojob=0;
else
    dojob=1;
end

if nargin<3
    njobs = 10;
end

if nargin<4
    nsecs = 2;
end

if dojob
    for ijob = 1:njobs
        fprintf('\n%s, %d/%d',jobname,ijob,njobs)
        jn = sprintf('@%s',jobname);
        mgsub({},@jobname,input,'qsub_optqsub_opts','-l h_vmem=16G')
        pause(nsecs)
    end
end
fprintf('\n')

%% sensorspace 
cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';

mgsub({},@fp_test_sensorspace_coh_allchans_group,{[],DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')

% for i_afs = 1:11
%     mgsub({},@fp_test_sensorspace_coh_allchans,{[],DIROUT,DIRLOG},'qsub_optqsub_opts','-l h_vmem=16G')
% end

%% touch lognames by hand 
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors
for id = 7
    for ichunk = 6
        logname = sprintf('%s_%d',patientID{id},ichunk);
        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
    end 
end


%% surrogates

cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
DIRLOG = '/home/bbci/data/haufe/Franziska/log/surrogate/';

for i=1:550
    mgsub({},@fp_surrogate_coh,{[],DIROUT,DIRLOG},'qsub_optqsub_opts','-l h_vmem=16G')
    pause(1)
end
%% thresholds (only ss)

DIROUT = '/home/bbci/data/haufe/Franziska/data/';
afs = {'theta','alpha','beta','gamma_low','gamma_high'};
abs_imag = {'abs','imag'};

cd ~/matlab/fp/

mgsub({},@fp_get_thresholds_ss_freq,{[],'abs',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')
mgsub({},@fp_get_thresholds_ss_freq,{[],'imag',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')
pause(2)

for i_afs = 1:numel(afs)
    for i_abs_imag = 1:numel(abs_imag)       
      mgsub({},@fp_get_thresholds_ss,{[],afs{i_afs},abs_imag{i_abs_imag},DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')  
      pause(2)
    end 
end 

%% clusters ss
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
afs = {'theta','alpha','beta','gamma_low','gamma_high'};
abs_imag = {'abs','imag'};
cd ~/matlab/fp/

for i_afs = 1:numel(afs)
    for i_abs_imag = 1:numel(abs_imag)       
      mgsub({},@fp_cluster_ss,{[],afs{i_afs},2,abs_imag{i_abs_imag},DIROUT},'qsub_optqsub_opts','-l h_vmem=16GB')  
      pause(2)
    end 
end 

mgsub({},@fp_cluster_ss_freq,{[],2,'abs',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')
mgsub({},@fp_cluster_ss_freq,{[],2,'imag',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')
pause(2)
mgsub({},@fp_cluster_ss_c_freq,{[],'abs',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G') %components fun 
mgsub({},@fp_cluster_ss_c_freq,{[],'imag',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')

%% cluster group
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
afs = {'theta','alpha','beta','gamma_low','gamma_high'};
abs_imag = {'abs','imag'};
cd ~/matlab/fp/

for i_afs = 1:numel(afs)
    for i_abs_imag = 1:numel(abs_imag)       
      mgsub({},@fp_cluster_g,{2, afs{i_afs},abs_imag{i_abs_imag},DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')  
      pause(2)
    end 
end 

mgsub({},@fp_cluster_g_freq,{2,'abs',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G') 
mgsub({},@fp_cluster_g_freq,{2,'imag',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')
pause(2)
mgsub({},@fp_cluster_g_c_freq,{'abs',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G') %components fun 
mgsub({},@fp_cluster_g_c_freq,{'imag',DIROUT},'qsub_optqsub_opts','-l h_vmem=16G')


%% eloreta filters

cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
DIRLOG = '/home/bbci/data/haufe/Franziska/log/filters_e/';

for i=1:11
    mgsub({},@fp_savefilter_eloreta,{[],DIROUT,DIRLOG},'qsub_optqsub_opts','-l h_vmem=16G')
    pause(2)
end






