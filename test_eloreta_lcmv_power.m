clear all

rng(5)
DIRFIG = '/Users/franziskapellegrini/Dropbox/Franziska/Data_MEG_Project/figures/mimsim_ana/mim_sim5/eLORETAvsLCMV6/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

DIROUT1=[];
params.iReg=1; %number of interacting voxels in interacting regions 
params.iInt = 1; %number of interacting regions
params.ilag = 2; %lag size
params.isnr = 0.9; %SNR
params.iss = 0.5; %noise mix
params.ip=[]; %paramenter configuration 
t=[]; %run time
params.ifilt='l';
params.dimred='p';

% sensor signal 
D = fp_get_Desikan(params.iReg);

%signal generation
fprintf('Signal generation... \n')
[sig,brain_noise,sensor_noise,L,iroi_seed, iroi_tar,D, fres, n_trials,filt,signal_sources] = ...
    fp_u(params,D,DIROUT1);
%%
noise = params.iss*brain_noise + (1-params.iss)*sensor_noise;   
noise_f = (filtfilt(filt.bband, filt.aband, noise'))';
noise = noise ./ norm(noise_f, 'fro');
%combine signal and noise
signal_sensor1 = params.isnr*sig + (1-params.isnr)*noise;
signal_sensor_f = (filtfilt(filt.bband, filt.aband, signal_sensor1'))';
signal_sensor1 = signal_sensor1 ./ norm(signal_sensor_f, 'fro');

%high-pass signal
signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
signal_sensor = signal_sensor / norm(signal_sensor, 'fro');

%reshape
signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);

% eloreta 
[n_sensors, l_epoch, n_trials] = size(signal_sensor);
L_backward = L(:, D.ind_cortex, :);
ndim = size(L_backward,3);
%%
reg_param = fp_eloreta_crossval(signal_sensor,L_backward,5);
A = squeeze(mkfilt_eloreta_v2(L_backward,reg_param));
A = permute(A,[1, 3, 2]);

%
source=[];
for idim = 1:3
    source(:,:,idim) = (squeeze(A(:,idim,:))'*signal_sensor(:,:))'.*10^9;
end 


% filtering 

low = [8 12];
% filters for band and highpass
[bband_low, aband_low] = butter(5, low/(fres/2));
source = filtfilt(bband_low, aband_low, source(:,:)); 

%

for is = 1:size(L_backward,2)
    va_ = var(squeeze(source(:, is,:)));
    sp(is) = sum(va_)';
end

%
load cm17;
load('processed_bs_wzb_90_2000/bs_results.mat')

data = zeros(2002,1);
data(D.ind_cortex)=sp;
% data = sp';
outname1 = [DIRFIG 'eLORETA_power_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data(in_normal_to_high), [0 max(abs(sp))], cm17a ,'Power', 0.3,outname1)


      


%% LCMV

cCS = cov(signal_sensor(:,:)');
reg = 0.05*trace(cCS)/length(cCS);
Cr = cCS + reg*eye(size(cCS,1));

[~, A_l] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
A_l = permute(A_l,[1, 3, 2]);

%
source_l = []; 
for idim = 1:3
    source_l(:,:,idim) = (squeeze(A_l(:,idim,:))'*signal_sensor(:,:))'.*10^9;
end 

% filtering 

low = [8 12];
% filters for band and highpass
[bband_low, aband_low] = butter(5, low/(fres/2));
source_l = filtfilt(bband_low, aband_low, source_l(:,:));  

%

for is = 1:size(L_backward,2)
    va_ = var(squeeze(source_l(:, is,:)));
    sp_l(is) = sum(va_)';
end

%
load cm17;
load('processed_bs_wzb_90_2000/bs_results.mat')

data_l = zeros(2002,1);
data_l(D.ind_cortex)=sp_l;
% data_l = sp_l';
outname1 = [DIRFIG 'LCMV_power_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data_l(in_normal_to_high), [0 max(abs(sp_l))], cm17a ,'Power', 0.3,outname1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% eloreta

% Loop over different pipelines

ipip = 1;

fprintf(['Testing pipeline ' num2str(ipip) '\n'])
% PCA

clear npcs  variance_explained
signal_roi = [];
empty_rois =[];
active_rois = [];

%loop over regions
for aroi = 1:D.nroi
    
    clear A_ signal_source
    A_ = A(:, :,D.ind_roi_cortex{aroi},:);
    
    
    %number of voxels at the current roi
    nvoxroi(aroi) = size(A_,3);
    A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
    
    %project sensor signal to voxels at the current roi (aroi)
    signal_source = A2{aroi}' * signal_sensor(:,:);
    
    
    clear s_
    
    
    s_ = logical(ones(size(signal_source,1),1));
    
    
    
    %do PCA
    clear signal_roi_ S
    [signal_roi_,S,~] = svd(double(signal_source(s_,:))','econ');
    
    
    npcs(aroi) = ipip; %fixed number of pcs
    
    
    %bring signal_roi to the shape of npcs x l_epoch x n_trials
    signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
    
end



% calculate connectivity scores

% calculate indices
clear inds PCA_inds
[inds, PCA_inds] = fp_npcs2inds(npcs);

output = {'MIM','MIC'};
conn = data2sctrgcmim(signal_roi, fres,30, 0,0, [], inds, output,0);

% extract measures out of the conn struct
[MIM_, MIC_, GC_, DIFFGC_, iCOH_, aCOH_] = fp_unwrap_conn(conn,numel(npcs),filt,PCA_inds);

data_mim = sum(MIM_,2);
outname1 = [DIRFIG 'eLORETA_netMIM_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data_mim, [0 max(abs(data_mim))], cm17a ,'netMIM score', 0.3,outname1)
%
% LCMV

% Loop over different pipelines

errorpipeline = [];
ipip = 1;

clear npcs  variance_explained
signal_roi_l = [];
empty_rois =[];
active_rois = [];

%loop over regions
for aroi = 1:D.nroi
    
    clear A_ signal_sourc
    A_ = A_l(:, :,D.ind_roi_cortex{aroi},:);
    
    %number of voxels at the current roi
    nvoxroi(aroi) = size(A_,3);
    A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
    
    %project sensor signal to voxels at the current roi (aroi)
    signal_source = A2{aroi}' * signal_sensor(:,:);
    
    clear s_
    
    s_ = logical(ones(size(signal_source,1),1));
    
    
    %do PCA
    clear signal_roi_ S
    [signal_roi_,S,~] = svd(double(signal_source(s_,:))','econ');
    
    
    npcs(aroi) = ipip; %fixed number of pcs
    
    
    
    %bring signal_roi to the shape of npcs x l_epoch x n_trials
    signal_roi_l = cat(1,signal_roi_l,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
    
    
end



% calculate connectivity scores

if ipip ~= 10 && ipip ~= 11 && ipip ~= 12
    
    if strcmp(params.ifilt(1),'c')
        %no calculations for empty rois
        npcs(empty_rois) = [];
    end
    
    % calculate indices
    clear inds PCA_inds
    [inds, PCA_inds] = fp_npcs2inds(npcs);
    
    % calculate MIM, MIC, TRGC and COH
    %             if ipip > 20
    %                 output = {'MIM','MIC','COH'};
    %             else
    output = {'MIM','MIC'};
    %             end
    conn = data2sctrgcmim(signal_roi_l, fres,30, 0,0, [], inds, output,0);
    
    % extract measures out of the conn struct
    [MIM_l, MIC_l, GC_, DIFFGC_, iCOH_, aCOH_] = fp_unwrap_conn(conn,numel(npcs),filt,PCA_inds);
    
    
end

%
close all
data_mim_l = sum(MIM_l,2);
outname1 = [DIRFIG 'LCMV_netMIM_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data_mim_l, [0 max(abs(data_mim_l))], cm17a ,'netMIM score', 0.3,outname1)
close all

[pr_mim_e(ipip)] = fp_pr(MIM_,iroi_seed,iroi_tar,1)
[pr_mim_l(ipip)] = fp_pr(MIM_l,iroi_seed,iroi_tar,1)


%

data_mim = squeeze(MIM_(iroi_seed,:));
outname1 = [DIRFIG 'eLORETA_seedMIM_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data_mim, [0 max(abs(data_mim))], cm17a ,'MIM score', 0.3,outname1)
close all

data_mim_l = squeeze(MIM_l(iroi_seed,:));
outname1 = [DIRFIG 'LCMV_seedMIM_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data_mim_l, [0 max(abs(data_mim_l))], cm17a ,'MIM score', 0.3,outname1)
close all

data_mim = squeeze(MIM_(iroi_tar,:));
outname1 = [DIRFIG 'eLORETA_tarMIM_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data_mim, [0 max(abs(data_mim))], cm17a ,'MIM score', 0.3,outname1)
close all

data_mim_l = squeeze(MIM_l(iroi_tar,:));
outname1 = [DIRFIG 'LCMV_tarMIM_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data_mim_l, [0 max(abs(data_mim_l))], cm17a ,'MIM score', 0.3,outname1)
close all
% corrs 


e = signal_roi(:,:);
l = signal_roi_l(:,:);
s = signal_sources';


for ii = 1:68
    
    for sr = 1:2 
        ce(ii,sr) = corr(e(ii,:)',s(sr,:)');
        cl(ii,sr) = corr(l(ii,:)',s(sr,:)');
    end
    
end


data = abs(squeeze(ce(:,1)));
outname1 = [DIRFIG 'eLORETA_seedcorr_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data, [0 max(abs(data))], cm17a ,'corr coeff', 0.3,outname1)
close all

data = abs(squeeze(ce(:,2)));
outname1 = [DIRFIG 'eLORETA_tarcorr_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data, [0 max(abs(data))], cm17a ,'.', 0.3,outname1)

data = abs(squeeze(cl(:,1)));
outname1 = [DIRFIG 'LCMV_seedcorr_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data, [0 max(abs(data))], cm17a ,'.', 0.3,outname1)
close all

data = abs(squeeze(cl(:,2)));
outname1 = [DIRFIG 'LCMV_tarcorr_snr0' num2str(params.isnr*10)];
allplots_cortex_BS(cortex_highres,data, [0 max(abs(data))], cm17a ,'.', 0.3,outname1)
close all