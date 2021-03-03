function [sig,brain_noise,sensor_noise, gt,L_save,iroi_seed,iroi_tar,D, fres, n_trials,filt] = fp_generate_signal_with_timestructure...
    (params,D,DIROUT1)

%if second condition of lag, then load parameters of first condition 
flag = true; 
if params.ip==6 
    params_save = params; 
    load(sprintf('%smim_lag/%d.mat',DIROUT1,params.iit));
    params = params_save; 
    clear params_save
    flag = false; 
end

if flag
    iroi_seed = randi(D.nroi,params.iInt,1);
    iroi_tar = randi(D.nroi,params.iInt,1);
    for ii = 1:params.iInt
        while any(iroi_seed==iroi_tar(ii))
            iroi_tar(ii) = randi(D.nroi,1,1);
        end
    end
end

%set parameters
fs = 100; % sampling rate
fres = fs; % number of frequency bins (= fres + 1)
Nmin = 3; % length of recording in minutes
N = Nmin*60*fs; % total number of samples
Lepo = 2*fres; % epoch length, should be consistent with fres
n_trials = N/Lepo; % number of epochs
% frqs = sfreqs(fres, fs); % freqs in Hz
iband = [8 12]; % frequency band of interaction in Hz
coupling_snr = 0.6; % coupling strength = SNR in interacting frequency band 
% ar_order = fres/5; % AR model order for TRGC estimation
%nboot = 30; % number of bootstrap iterations

% indices of interacting frequencies
% band_inds = find(frqs >= iband(1) & frqs <= iband(2));

% filters for band and highpass
[bband, aband] = butter(2, iband/fs*2);
[bhigh, ahigh] = butter(2, 1/fs*2, 'high');
filt.aband = aband;
filt.bband = bband;
filt.iband = iband;

%set random small or large lag
if params.ilag == 1
    lag = randi([0, 5],params.iInt*params.iReg,1);
else
    lag = randi([5, 20],params.iInt*params.iReg,1);
end
% lag_ms = fs/lag; % lag in ms

%% indices of signal and noise 
sig_ind = [];
for ii = 1:params.iReg
    sig_ind = [sig_ind; (iroi_seed.*params.iReg)-(ii-1), (iroi_tar.*params.iReg)-(ii-1)];
end

noise_ind = setdiff(1:params.iReg*D.nroi,sig_ind(:));

%% generate interacting sources 

if flag
    s1 = randn(N, params.iReg*params.iInt);
    s1 = filtfilt(bband, aband, s1);  
end

%save this state of s1 for later use
if params.ip==1
    s1_save = s1;
end

for ii = 1:params.iInt*params.iReg
    s2(:,ii) = circshift(squeeze(s1(:,ii)), lag(ii));
end

s1 = cat(2,s1,s2);
s1 = s1 / norm(s1, 'fro');

% pink background noise is added
if flag
    backg = mkpinknoise(N, params.iInt*params.iReg*2, 1);
    backgf = filtfilt(bband, aband, backg);
    % normalization is done w.r.t. interacting band
    backg = backg / norm(backgf, 'fro');
end

signal_sources = coupling_snr*s1 + (1-coupling_snr)*backg;

%% non-interacting sources

if flag
    noise_sources = mkpinknoise(N, params.iReg*D.nroi-(params.iReg*params.iInt*2), 1);
end

%% combine interacting and non-interacting sources and generate ground truth imaginary coherence

% signal_gt = zeros(N,params.iReg*D.nroi);
% 
% signal_gt(:,sig_ind(:,1)) = signal_sources(:,1:params.iReg*params.iInt); %seed voxels
% signal_gt(:,sig_ind(:,2)) = signal_sources(:,params.iReg*params.iInt+1:end);%target voxels 
% 
% noise_sources_f = filtfilt(bband, aband, noise_sources);
% noise_sources1 = noise_sources / norm(noise_sources_f, 'fro');
% signal_gt(:,noise_ind) = noise_sources1;
% 
% signal_gt = permute(signal_gt,[2 1]);
% 
% % generate ground truth imaginary coherence
% signal_gt_trials = reshape(signal_gt, params.iReg*D.nroi, Lepo, n_trials);
% CS_gt = tsdata_to_cpsd_fast(signal_gt_trials,fres,'WELCH');
% CS_gt(:,:,[1])=[];
% for ifreq = 1: fres
%     clear pow
%     pow = real(diag(CS_gt(:,:,ifreq)));
%     gt1(:,:,ifreq) = CS_gt(:,:,ifreq)./ sqrt(pow*pow');
% end
% if params.iReg~=1 %mim across two voxels of one region 
%     [gt_mic,gt_mim]= fp_mim(gt1,repmat(params.iReg,D.nroi,1));
%     gt.mic = gt_mic; 
%     gt.mim = gt_mim; 
% else 
%     gt.mic = abs(imag(gt1));
%     gt.mim = abs(imag(gt1)); 
% end
% clear gt1 gt_mic gt_mim

%% leadfield for forward model

L_save = D.leadfield;
L3 = L_save(:, D.sub_ind_cortex, :);
normals = D.normals(D.sub_ind_cortex,:)'; 
for is = 1:numel(D.sub_ind_cortex)
    L_mix(:,is) = squeeze(L3(:,is,:))*squeeze(normals(:,is));
end

%select signal L and noise L 
L_sig = L_mix(:,sig_ind);
L_noise = L_mix(:,noise_ind);

%% project to sensors and add white noise 

%signal
sig = L_mix(:,sig_ind) * signal_sources';
sig_f = (filtfilt(bband, aband, sig'))';
sig = sig ./ norm(sig_f, 'fro'); 

%brain noise
if flag
    brain_noise = L_mix(:,noise_ind) * noise_sources';
    brain_noise_f = (filtfilt(bband, aband, brain_noise'))';
    brain_noise = brain_noise ./ norm(brain_noise_f, 'fro');
end

%white noise
if flag
    sensor_noise = randn(size(sig));
    sensor_noise_f = (filtfilt(bband, aband, sensor_noise'))';
    sensor_noise = sensor_noise ./ norm(sensor_noise_f, 'fro');
end

if params.ip==1
    fprintf('Saving lag stuff... \n')
    dir1 =  sprintf('%smim_lag/',DIROUT1);
    if ~exist(dir1); mkdir(dir1); end
    outname = sprintf('%smim_lag/%d.mat',DIROUT1,params.iit);
    s1 = s1_save;
    save(outname,'iroi_seed','iroi_tar','D','sensor_noise','brain_noise','backg','s1','-v7.3')   
end
