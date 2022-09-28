function fp_eval_lagsim
%Simulation on the comparison between different interaction delays.

fp_addpath

DIRLOG ='/home/bbci/data/haufe/Franziska/log/mim_sim5_lag/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

rng('shuffle')

%%
%prevent array jobs to start at exactly the same time
iit = str2num(getenv('SGE_TASK_ID'))

clear nInteractions nRegionInts SNR noise_mix nlag filtertype hemisym

ip = 1;
[nInteractions,nRegionInts,SNR,noise_mix,~,filtertype,~] = fp_get_params(ip);

iInt = nInteractions;
iReg = nRegionInts;
isnr = SNR;
iss = noise_mix;
ifilt = filtertype{1};

%logname
logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_filt%s_iter%d'...
    ,iInt,iReg,isnr*10,iss*10,ifilt,iit);

%create logfile for parallelization
if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
    eval(sprintf('!touch %s%s_work',DIRLOG,logname))
    fprintf('Working on %s. \n',logname)
    
    params.iInt = iInt;
    params.iReg = iReg;
    params.isnr = isnr;
    params.iss = iss;
    params.ifilt = ifilt;
    params.iit = iit;
    params.ip = ip;
    params.logname = logname;
    
    fp_lag_sim(params)
    
    eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
    
end  %work_done



