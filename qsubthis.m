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
        mgsub({},@jobname,input,'qsub_opts','-l h_vmem=16G')
        pause(nsecs)
    end
end
fprintf('\n')

%% sensorspace 
cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
afs = {'theta','alpha','beta','gamma_low','gamma_high'};

for i_afs= 1:4
    mgsub({},@fp_test_sensorspace_coh_allchans_group,{[],afs{i_afs},DIROUT},'qsub_opts','-l h_vmem=16G')
end

% for i_afs = 1:11
%     mgsub({},@fp_test_sensorspace_coh_allchans,{[],DIROUT,DIRLOG},'qsub_opts','-l h_vmem=16G')
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

for i=1:10
    
    mgsub({},@fp_surrogate_coh,{[],DIROUT,DIRLOG},'qsub_opts','-l h_vmem=16G')
    pause(1)
end

%% surrogates eloreta 

cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
DIRLOG = '/home/bbci/data/haufe/Franziska/log/surrogate/eloreta/';

for i=1:500
    
    mgsub({},@fp_surrogate_coh_e,{[],DIROUT,DIRLOG},'qsub_opts','-l h_vmem=16G')
    pause(1)
end

%% surrogates eloreta 2D

cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
DIRLOG = '/home/bbci/data/haufe/Franziska/log/surrogate/eloreta/2D/';

for i=1:40
    
    mgsub({},@fp_surrogate_coh_e2D,{[],DIROUT,DIRLOG},'qsub_opts','-l h_vmem=16G')
    pause(1)
end
%% thresholds (only ss)

DIROUT = '/home/bbci/data/haufe/Franziska/data/';
afs = {'theta','alpha','beta','gamma_low','gamma_high'};
abs_imag = {'abs','imag'};

cd ~/matlab/fp/

mgsub({},@fp_get_thresholds_ss_freq,{[],'abs',DIROUT},'qsub_opts','-l h_vmem=16G')
mgsub({},@fp_get_thresholds_ss_freq,{[],'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
pause(1)

for i_afs = 1:numel(afs)
    for i_abs_imag = 1:numel(abs_imag)       
      mgsub({},@fp_get_thresholds_ss,{[],afs{i_afs},abs_imag{i_abs_imag},DIROUlllT},'qsub_opts','-l h_vmem=16G')  
      pause(1)
    end 
end 

%% clusters ss
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
afs = {'theta','alpha','beta','gamma_low','gamma_high'};
abs_imag = {'abs','imag'};
cd ~/matlab/fp/

% for i_afs = 1:numel(afs)
%     for i_abs_imag = 1:numel(abs_imag)       
%       mgsub({},@fp_cluster_ss,{[],afs{i_afs},2,abs_imag{i_abs_imag},DIROUT},'qsub_opts','-l h_vmem=16GB')  
%       pause(2)
%     end 
% end 
% 
% mgsub({},@fp_cluster_ss_freq,{[],2,'abs',DIROUT},'qsub_opts','-l h_vmem=16G')
% mgsub({},@fp_cluster_ss_freq,{[],2,'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
% pause(2)
% mgsub({},@fp_cluster_ss_c_freq,{[],'abs',DIROUT},'qsub_opts','-l h_vmem=16G') %components fun 
% mgsub({},@fp_cluster_ss_c_freq,{[],'imag',DIROUT},'qsub_opts','-l h_vmem=16G')

%% cluster group
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
afs = {'theta','alpha','beta','gamma_low','gamma_high'};
abs_imag = {'abs','imag'};
cd ~/matlab/fp/
% 
% for i_afs = 1:numel(afs)
%     for i_abs_imag = 1:numel(abs_imag)       
%       mgsub({},@fp_cluster_g,{2, afs{i_afs},abs_imag{i_abs_imag},DIROUT},'qsub_opts','-l h_vmem=16G')  
%       pause(2)
%     end 
% end 
% 
mgsub({},@fp_cluster_g_freq,{2,'abs',DIROUT},'qsub_opts','-l h_vmem=16G') 
mgsub({},@fp_cluster_g_freq,{2,'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
pause(2)
mgsub({},@fp_cluster_g_c_freq,{'abs',DIROUT},'qsub_opts','-l h_vmem=16G') %components fun 
mgsub({},@fp_cluster_g_c_freq,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')

%only julian's subs
for i_afs = 1:numel(afs)
    for i_abs_imag = 1:numel(abs_imag)       
      mgsub({},@fp_cluster_g_j2,{2, afs{i_afs},abs_imag{i_abs_imag},DIROUT},'qsub_opts','-l h_vmem=16G')  
      pause(2)
    end 
end 

%only julian's subs
mgsub({},@fp_cluster_g_freq_j,{2,'abs',DIROUT},'qsub_opts','-l h_vmem=16G') 
mgsub({},@fp_cluster_g_freq_j,{2,'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
pause(2)
mgsub({},@fp_cluster_g_c_freq_j,{'abs',DIROUT},'qsub_opts','-l h_vmem=16G') %components fun 
mgsub({},@fp_cluster_g_c_freq_j,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')


%new cluster method
% mgsub({},@fp_cluster_g_freq_test,{2,'abs',DIROUT},'qsub_opts','-l h_vmem=16G') 
mgsub({},@fp_cluster_g_freq_test,{2,'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
pause(1)
% mgsub({},@fp_cluster_g_c_freq_test,{'abs',DIROUT},'qsub_opts','-l h_vmem=16G') %components fun 
mgsub({},@fp_cluster_g_c_freq_test,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')


%new cluster method + 2D
% mgsub({},@fp_cluster_g_freq_test,{2,'abs',DIROUT},'qsub_opts','-l h_vmem=16G') 
mgsub({},@fp_cluster_g_freq_test2D,{2,'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
pause(1)
% mgsub({},@fp_cluster_g_c_freq_test,{'abs',DIROUT},'qsub_opts','-l h_vmem=16G') %components fun 
mgsub({},@fp_cluster_g_c_freq_test2D,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')

%more parameters varied 
mgsub({},@fp_cluster_g_c_freq_test_0005,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
mgsub({},@fp_cluster_g_c_freq_test_fwf,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')


%dics 
mgsub({},@fp_cluster_g_freq_test_d,{2,'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
mgsub({},@fp_cluster_g_c_freq_test_d,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
mgsub({},@fp_cluster_g_c_freq_test_0005_d,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')
mgsub({},@fp_cluster_g_c_freq_test_fwf_d,{'imag',DIROUT},'qsub_opts','-l h_vmem=16G')


%% cluster proto function 
cd ~/matlab/fp/

DIROUT = '/home/bbci/data/haufe/Franziska/data/';
clusterfun = {'c','f'};
abs_imag = {'abs','imag'};
testmethod = {'s','t'};
alpha = [0.005, 0.001, 0.0005];
fwf = [0, 1];
j = [0, 1];
filtertype = {'e','e2D','d'};

%specific job 
mgsub({},@fp_cluster_proto,{DIROUT,filtertype{3}, clusterfun{2},abs_imag{2},testmethod{1},alpha(2)...
                            ,fwf(1), j(2)},'qsub_opts','-l h_vmem=16G')

%all
for ft = 1:3
    for c = 1:2
        for ai = 1:2
            for t = 1:2
                for al = 1:3
                    for f = 1:2
                        for jn = 1:2
                            mgsub({},@fp_cluster_proto,{DIROUT,filtertype{ft}, clusterfun{c},abs_imag{ai},testmethod{t},alpha(al)...
                                ,fwf(f), j(jn)},'qsub_opts','-l h_vmem=16G')
                            pause(1)
                        end
                    end
                end
            end 
        end 
    end 
end

%% eloreta filters

cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
DIRLOG = '/home/bbci/data/haufe/Franziska/log/filters_e/';

for i=1:11
    mgsub({},@fp_savefilter_eloreta,{[],DIROUT,DIRLOG},'qsub_opts','-l h_vmem=16G')
    pause(2)
end

%% filters
cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';

DIRLOG = '/home/bbci/data/haufe/Franziska/log/filters_dics/';

for i=1:11
    mgsub({},@fp_savefilter_temp,{[],DIROUT,DIRLOG},'qsub_opts','-l h_vmem=16G')
    pause(2)
end

%% gc_pipeline

cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';

mgsub({},@fp_gc_pipeline,{[],DIROUT},'qsub_opts','-l h_vmem=16G')

%% megmeg_pipeline

cd ~/matlab/fp/
DIROUT = '/home/bbci/data/haufe/Franziska/data/';

mgsub({},@fp_megmeg_pipeline,{[],DIROUT},'qsub_opts','-l h_vmem=16G')

%% mim struct sim 
clear all
% fp_addpath
cd ~/matlab/fp/
nit = 100;
varyParam = [1:7];%defaults, snr, interactions, everything else 

for ip =7
    mgsub({},@fp_eval_mim_struct_sim,{ip},'qsub_opts',['-l h_vmem=16G -t 1-' num2str(nit) ]) %' -q 2jobs'
    pause(2)
end

%%

clear all
% fp_addpath
cd ~/matlab/fp/

mgsub({},@fp_test,{},'qsub_opts',['-l h_vmem=32G -t 1-3'])
