


%default paramenters
params.iInt = 1;
params.iReg = 1;
params.isnr = 0.5;
params.iss = 0.5;
params.ilag = 2;
params.ifilt = 'l';
params.ihemi = 0;
params.iit = 1;
params.ip = 1;

logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d_iter%d'...
    ,params.iInt,params.iReg,params.isnr*10,params.iss*10, params.ilag,params.ifilt,params.ihemi,params.iit);
params.logname = logname;

%% Option 1: call the whole simulation from scratch:

%calculates all pipelines (takes some time)
fp_mim_struct_sim(params)

%% Option2: Load existing results file and re-calculate the pipelines of interest 

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim/';
inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d_iter%d'...
        ,params.iInt,params.iReg,params.isnr*10,params.iss*10, params.ilag,params.ifilt,params.ihemi,params.iit);
    
load([DIRIN inname '.mat'])
clear mic mim 

%%
% mode1 is either a number for fixed pcs (an integer), or 'max' (select npcs = rank of
%region data), or 'percent' (select npcs that 90% of the variance is
%preserved), or 'case2' (mim only to pool dimensions, then summation), or
%'baseline', or 'all'
mode1='baseline';
[mic1, mim1, u1] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,mode1);

%In fp_get_mim, also the npcs are calculated and the PCA is done. 
%The selection of PCs and the calculation of the mim is then done in the 
%subfunction fp_compute_mode_mim. 

hist(mic1(:))
figure; 
hist(mim1(:))

%% 

[mic2, mim2, u2] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,mode1);
hist(mic2(:))
figure
hist(mim2(:))

%%
norm(mic1(:)-mic2(:),'fro')
norm(mim1(:)-mim2(:),'fro')

