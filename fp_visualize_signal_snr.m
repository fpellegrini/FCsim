
% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

rng(5)

params.iReg=1; %number of interacting voxels in interacting regions
params.iInt = 1; %number of interacting regions
params.ilag = 2; %lag size
params.isnr1 = 0.3; %SNR
params.isnr2 = 0.6; %SNR
params.isnr3 = 0.9; %SNR
params.iss = 0.5; %noise mix
params.ip=3; %paramenter configuration
t=[]; %run time
params.ifilt='l';
params.dimred='p';

% sensor signal
D = fp_get_Desikan(params.iReg);
no_reload = true;

%% set parameters

fs = 100; % sampling rate
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
    iroi_seed = 11;%randperm(D.nroi,params.iInt)';
    iroi_tar = 49;%randperm(D.nroi,params.iInt)';
    
    %be sure that no region is selected twice
    for ii = 1:params.iInt
        while any(iroi_seed==iroi_tar(ii)) || sum(iroi_tar==iroi_tar(ii))>1
            iroi_tar(ii) = randi(D.nroi,1,1);
        end
    end
end

%set random small or large lag
if params.ilag == 1
    lag = randi([1, 5],params.iInt*params.iReg,1);
else
    lag = randi([6, 20],params.iInt*params.iReg,1);
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
    s2(:,ii) = circshift(squeeze(s1(:,ii)), lag(ii));
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
L=L_save;
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
    try
        brain_noise = L_noise * noise_sources';
        brain_noise_f = (filtfilt(bband, aband, brain_noise'))';
        brain_noise = brain_noise ./ norm(brain_noise_f, 'fro');
    catch
        iroi_seed
        iroi_tar
    end
    
end

%white noise on sensor level (sensor noise)
if no_reload
    sensor_noise = randn(size(sig));
    sensor_noise_f = (filtfilt(bband, aband, sensor_noise'))';
    sensor_noise = sensor_noise ./ norm(sensor_noise_f, 'fro');
end

%combine noise sources
noise = params.iss*brain_noise + (1-params.iss)*sensor_noise;
noise_f = (filtfilt(filt.bband, filt.aband, noise'))';
noise = noise ./ norm(noise_f, 'fro');

%%
%
clear signal_sensor1 signal_sensor_f signal_sensor
%combine signal and noise
signal_sensor1 = params.isnr1*sig + (1-params.isnr1)*noise;
signal_sensor_f = (filtfilt(filt.bband, filt.aband, signal_sensor1'))';
signal_sensor1 = signal_sensor1 ./ norm(signal_sensor_f, 'fro');

%high-pass signal
signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
signal_sensor = signal_sensor / norm(signal_sensor, 'fro');

%reshape
signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
[n_sensors, l_epoch, n_trials] = size(signal_sensor);

%
ya = [-70 -49];
u=[];
CS = tsdata_to_cpsd_fast(signal_sensor,fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);


figsiz = [7,30];
figure
figone(figsiz(1),figsiz(2))
subplot(1,3,1)
plot(u')
xlim([2 100])
ylim([-80 -45])
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)
title('-7.4 dB')

xlabel('Frequency (Hz)','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on

%%%%
clear signal_sensor1 signal_sensor_f signal_sensor
%combine signal and noise
signal_sensor1 = params.isnr2*sig + (1-params.isnr2)*noise;
signal_sensor_f = (filtfilt(filt.bband, filt.aband, signal_sensor1'))';
signal_sensor1 = signal_sensor1 ./ norm(signal_sensor_f, 'fro');

%high-pass signal
signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
signal_sensor = signal_sensor / norm(signal_sensor, 'fro');

%reshape
signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
[n_sensors, l_epoch, n_trials] = size(signal_sensor);

%

u=[];
CS = tsdata_to_cpsd_fast(signal_sensor,fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

subplot(1,3,2)
plot(u')
xlim([2 100])
ylim([-80 -45])
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)

xlabel('Frequency (Hz)','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on
title('3.5 dB')

%%%%
clear signal_sensor1 signal_sensor_f signal_sensor
%combine signal and noise
signal_sensor1 = params.isnr3*sig + (1-params.isnr3)*noise;
signal_sensor_f = (filtfilt(filt.bband, filt.aband, signal_sensor1'))';
signal_sensor1 = signal_sensor1 ./ norm(signal_sensor_f, 'fro');

%high-pass signal
signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
signal_sensor = signal_sensor / norm(signal_sensor, 'fro');

%reshape
signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
[n_sensors, l_epoch, n_trials] = size(signal_sensor);

%
figsiz = [7,20];
ya = [-70 -49];
u=[];
CS = tsdata_to_cpsd_fast(signal_sensor,fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

subplot(1,3,3)
plot(u')
xlim([2 100])
ylim([-80 -45])
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)

xlabel('Frequency (Hz)','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on
title('19.1 dB')

%%

saveas(gca,'~/Desktop/signal_sensor_snr.png','png')
print('~/Desktop/signal_sensor_snr.eps','-depsc')
close all
