function [sig,brain_noise,sensor_noise, L_save,iroi_seed,iroi_tar,D, fres, n_trials,filt] = fp_generate_signal_lagsim...
    (params,D,DIROUT1,lag)

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

%%
no_reload = true;
if lag > 1 && params.ip~=1
    %reload data from ip1 to keep them constant and only vary the lag size
    params_save = params;
    load(sprintf('%smim_lag/%d.mat',DIROUT1,params.iit));
    params = params_save;
    clear params_save
    no_reload = false;
end

%% set parameters

fs = 500; % sampling rate
fres = fs; % number of frequency bins (= fres + 1)
Nmin = 3; % length of recording in minutes
N = Nmin*60*fs; % total number of samples
Lepo = 2*fres; % epoch length, should be consistent with fres
n_trials = N/Lepo; % number of epochs
frqs = sfreqs(fres, fs); % freqs in Hz
iband = [8 12]; % frequency band of interaction in Hz
coupling_snr = 0.6; % coupling strength = SNR in interacting frequency band 
band_inds = find(frqs >= iband(1) & frqs <= iband(2)); % indices of interacting frequencies

% filters for band and highpass
[bband, aband] = butter(2, iband/fs*2);
[bhigh, ahigh] = butter(2, 1/fs*2, 'high');
filt.aband = aband;
filt.bband = bband;
filt.ahigh = ahigh;
filt.bhigh = bhigh;
filt.band_inds = band_inds;

if no_reload    
    %set seed and target regions 
    iroi_seed = randperm(D.nroi,params.iInt)';
    iroi_tar = randperm(D.nroi,params.iInt)';
    
    %be sure that no region is selected twice 
    for ii = 1:params.iInt
        while any(iroi_seed==iroi_tar(ii))
            iroi_tar(ii) = randi(D.nroi,1,1);
        end
    end   
end

%% indices of signal and noise 

sig_ind = [];
for ii = 1:params.iReg
    sig_ind = [sig_ind; (iroi_seed.*params.iReg)-(ii-1), (iroi_tar.*params.iReg)-(ii-1)];
end

noise_ind = setdiff(1:params.iReg*D.nroi,sig_ind(:));

%% generate interacting sources 

if no_reload
    %generate filtered white noise at seed voxels 
    s1 = randn(N, params.iReg*params.iInt);
    s1 = filtfilt(bband, aband, s1);  
end

% if ip1, save this state of s1 for ip6
if params.ip==1
    s1_save = s1;
end

for ii = 1:params.iInt*params.iReg
    %activity at target voxels is a shifted version of the seed voxels
    s2(:,ii) = circshift(squeeze(s1(:,ii)), lag);
end

%concenate seed and target voxel activity
s1 = cat(2,s1,s2);
s1 = s1 / norm(s1, 'fro');

% pink background noise is added
if no_reload
    backg = mkpinknoise(N, params.iInt*params.iReg*2, 1);
    backgf = filtfilt(bband, aband, backg);
    % normalization is done w.r.t. interacting band
    backg = backg / norm(backgf, 'fro');
end

%combine signal and background noise 
signal_sources = coupling_snr*s1 + (1-coupling_snr)*backg;

%% non-interacting sources

if no_reload
    %activity at all voxels but the seed and target voxels 
    noise_sources = mkpinknoise(N, params.iReg*D.nroi-(params.iReg*params.iInt*2), 1);
end

%% leadfield for forward model

L_save = D.leadfield;
L3 = L_save(:, D.sub_ind_cortex, :); % select only voxels that belong to a region 

% multiply with normal direction to get from three to one dipole dimension 
normals = D.normals(D.sub_ind_cortex,:)'; 
for is = 1:numel(D.sub_ind_cortex)
    L_mix(:,is) = squeeze(L3(:,is,:))*squeeze(normals(:,is));
end

%select signal L and noise L 
L_sig = L_mix(:,sig_ind);
L_noise = L_mix(:,noise_ind);

%% project to sensors and generate white noise 

%signal on sensor level 
sig = L_sig * signal_sources';
sig_f = (filtfilt(bband, aband, sig'))';
sig = sig ./ norm(sig_f, 'fro'); 

%brain noise on sensor level 
if no_reload
    brain_noise = L_noise * noise_sources';
    brain_noise_f = (filtfilt(bband, aband, brain_noise'))';
    brain_noise = brain_noise ./ norm(brain_noise_f, 'fro');
end

%white noise on sensor level (sensor noise) 
if no_reload
    sensor_noise = randn(size(sig));
    sensor_noise_f = (filtfilt(bband, aband, sensor_noise'))';
    sensor_noise = sensor_noise ./ norm(sensor_noise_f, 'fro');
end

% if lag 1, save this state
if params.ip==1 && lag==1
    fprintf('Saving lag stuff... \n')
    dir1 =  sprintf('%smim_lag/',DIROUT1);
    if ~exist(dir1); mkdir(dir1); end
    outname = sprintf('%smim_lag/%d.mat',DIROUT1,params.iit);
    s1 = s1_save;
    save(outname,'iroi_seed','iroi_tar','D','sensor_noise','brain_noise','backg','s1','-v7.3')   
end
