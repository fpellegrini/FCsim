function fp_figure1_and_5
%Plots two figures:
%Figure 1: Visualize different parts of the signal
%Figure 5: plot ground truth vs MIM scores vs p-values
%null distribution for p-values generated with fp_simulation_stats, and
%fp_submit_simulation_stats_nsg (which calls fp_simulation_stats2 on nsg
%cluster)

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

rng('default')
rng(5)

%input path for p-values 
DIRIN = '~/Dropbox/Franziska/Data_MEG_Project/simulation_stats/5/';

params.iReg=1; %number of interacting voxels in interacting regions
params.iInt = 1; %number of interacting regions
params.ilag = 2; %lag size
params.isnr = 0.6; %SNR
params.iss = 0.5; %noise mix
params.ip=3; %paramenter configuration
ipip=1;
params.ifilt='l';
params.dimred='p';
params.iit = 1;

% sensor signal
D = fp_get_Desikan(params.iReg);
no_reload = true;
ipip=3;

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

sources = zeros(N,D.nroi);
sources(:,sig_ind) = signal_sources;
sources(:,noise_ind) = noise_sources;

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
%combine signal and noise
signal_sensor1 = params.isnr*sig + (1-params.isnr)*noise;
signal_sensor_f = (filtfilt(filt.bband, filt.aband, signal_sensor1'))';
signal_sensor1 = signal_sensor1 ./ norm(signal_sensor_f, 'fro');

%high-pass signal
signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
signal_sensor = signal_sensor / norm(signal_sensor, 'fro');

%reshape
signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
[n_sensors, l_epoch, n_trials] = size(signal_sensor);


%select only voxels that belong to any roi
L_backward = L(:, D.ind_cortex, :);
ndim = size(L_backward,3); % 3 spatial dimensions

cCS = cov(signal_sensor(:,:)');
reg = 0.05*trace(cCS)/length(cCS);
Cr = cCS + reg*eye(size(cCS,1));

[~, A] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
A = permute(A,[1, 3, 2]);

%% dimensionality reduction

clear npcs  variance_explained
signal_roi = [];

%loop over regions
for aroi = 1:D.nroi
    
    clear A_ signal_source
    A_ = A(:, :,D.ind_roi_cortex{aroi},:);
    
    %number of voxels at the current region
    nvoxroi(aroi) = size(A_,3);
    A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
    
    %project sensor signal to voxels at the current roi (aroi)
    signal_source = A2{aroi}' * signal_sensor(:,:);
    
    %do PCA
    clear signal_roi_ S
    [signal_roi_,S,~] = svd(double(signal_source(:,:))','econ');
    
    % variance explained
    vx_ = cumsum(diag(S).^2)./sum(diag(S).^2);
    invx = 1:min(length(vx_), n_sensors);
    varex = vx_(invx);
    
    %fixed number of pcs
    npcs(aroi) = ipip;
   
    %bring signal_roi to the shape of npcs x l_epoch x n_trials
    signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));    
end

%% calculate MIM of original source activity for figure 5

% calculate indices
clear inds1 PCA_inds1
npcs1 = ones(1,D.nroi);
[inds1, PCA_inds1] = fp_npcs2inds(npcs1);

%true MIM 
output = {'MIM'};
sources1 = reshape(sources',D.nroi,size(signal_roi,2),size(signal_roi,3));
conn1 = data2sctrgcmim(sources1, fres, 30, 0,0, [], inds1, output);
[MIM_o, ~, ~, ~, ~,~] = fp_unwrap_conn(conn1,D.nroi,filt,PCA_inds1);

%% calculate MIM of reconstructed source activity for figure 5

% calculate indices
clear inds PCA_inds
[inds, PCA_inds] = fp_npcs2inds(npcs);

%true MIM 
output = {'MIM'};
conn = data2sctrgcmim(signal_roi, fres, 30, 0,0, [], inds, output);
[MIM_t, ~, ~, ~, ~,~] = fp_unwrap_conn(conn,D.nroi,filt,PCA_inds);

[pr] = fp_pr(MIM_t,iroi_seed,iroi_tar,1);  

%% calculate p-values for figure 5

load([DIRIN 'signal.mat'])

mim_s = [];
for ii = 1:100 
    try
    clear MIM_s
    load([DIRIN 'result_' num2str(ii) '.mat'])
    mim_s = cat(3,mim_s,MIM_s);
    end
end


for iroi = 1:D.nroi 
    for jroi = 1:D.nroi
        MIM_p(iroi,jroi) = sum(squeeze(mim_s(iroi,jroi,:))>MIM_t(iroi,jroi))/size(mim_s,3);
    end
end

ind_triu = logical(triu(ones(D.nroi,D.nroi)));
[p_fdr,~] = fdr(MIM_p(ind_triu),0.05);
MIM_p(MIM_p>p_fdr)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure 1 

%figure 1a
figsiz = [7,20];
figure
figone(figsiz(1),figsiz(2))

subplot(1,2,1)
plot(s1(20:220,:))
xlim([1 100])
xlabel('Time (msec)','FontSize',16)
ylabel('a.u.','FontSize',16)
grid on
xticks = 20:20:100;
xticklab = 200:200:1000;
set(gca,'xtick',xticks,'xticklabel',xticklab)

subplot(1,2,2)
u=[];
CS = tsdata_to_cpsd_fast(s1',fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end
u=10*log10(u);

plot(u')
xlim([2 100])
% ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)
xlabel('Frequency (Hz)','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on

saveas(gca,'~/Desktop/puresig.png','png')
print('~/Desktop/puresig.eps','-depsc')
close all


%figure 1b
figure
figone(figsiz(1),figsiz(2))
subplot(1,2,1)
plot(signal_sources(20:220,:))
xlim([1 100])
xlabel('Time (msec)','FontSize',16)
ylabel('a.u.','FontSize',16)
grid on
xticks = 20:20:100;
xticklab = 200:200:1000;
set(gca,'xtick',xticks,'xticklabel',xticklab)

ya = [-50 -20];
u=[];
CS = tsdata_to_cpsd_fast(signal_sources',fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end
u=10*log10(u);

subplot(1,2,2)
plot(u')
xlim([2 100])
ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)

xlabel('Frequency (Hz)','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on
saveas(gca,'~/Desktop/interactive.png','png')
print('~/Desktop/interaction.eps','-depsc')
close all


%figure 1c
figure
figone(figsiz(1),figsiz(2))

subplot(1,2,1)
plot(noise_sources(20:220,1))
xlim([1 100])
xlabel('Time (msec)','FontSize',16)
ylabel('a.u.','FontSize',16)
grid on
xticks = 20:20:100;
xticklab = 200:200:1000;
set(gca,'xtick',xticks,'xticklabel',xticklab)

ya = [-40 -15];
u=[];
CS = tsdata_to_cpsd_fast(noise_sources(:,1)',fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

subplot(1,2,2)
plot(u')
xlim([2 100])
ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)
xlabel('Frequency (Hz)','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on

saveas(gca,'~/Desktop/noninteractive.png','png')
print('~/Desktop/noninteractive.eps','-depsc')
close all


%figure 1d
ya = [-70 -49];
u=[];
CS = tsdata_to_cpsd_fast(signal_sensor,fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

figure
figone(figsiz(1),(figsiz(2)/2)-1)
plot(u')
xlim([2 100])
% ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)

xlabel('Frequency (Hz)','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on
saveas(gca,'~/Desktop/signal_sensor.png','png')
print('~/Desktop/signal_sensor.eps','-depsc')
close all

%figure 1e
%select first roi 
for iroi = 1:length(PCA_inds)   
    signal_roi_1stPC(iroi,:,:)= signal_roi(PCA_inds{iroi}(1),:,:);
end
ya = [-70 -49];
u=[];
CS = tsdata_to_cpsd_fast(signal_roi_1stPC(:,:),fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

figure
figone(figsiz(1),(figsiz(2)/2)-1)
plot(u')
xlim([2 100])
% ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)
xlabel('Frequency (Hz)','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on

saveas(gca,'~/Desktop/signal_roi.png','png')
print('~/Desktop/signal_roi.eps','-depsc')
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 5 

%figure 5a
figure 
figone(11,14)
imagesc(MIM_o)
set(gca,'FontSize',16)
caxis([0 0.15])
a=colorbar;
ylabel(a,'MIM score','FontSize',16,'Rotation',270);
a.FontSize = 16;
a.Label.Position(1) = 4.4;
xlabel('Region number','FontSize',16)
ylabel('Region number','FontSize',16)
set(gca,'FontSize',16)
saveas(gca,'~/Desktop/gt_connectome.png','png')
print('~/Desktop/gt_connectome.eps','-depsc')
close all


%figure 5b
figure 
figone(11,14)
imagesc(MIM_t)
set(gca,'FontSize',16)
caxis([0 0.15])
a=colorbar;
ylabel(a,'MIM score','FontSize',16,'Rotation',270);
a.FontSize = 16;
a.Label.Position(1) = 4.4;
xlabel('Region number','FontSize',16)
ylabel('Region number','FontSize',16)
set(gca,'FontSize',16)
saveas(gca,'~/Desktop/reconstructed_connectome.png','png')
print('~/Desktop/reconstructed_connectome.eps','-depsc')
close all


%figure 5c
figure 
figone(11,14)
imagesc(-log10(MIM_p))
set(gca,'FontSize',16)
caxis([0 4])
a=colorbar;
ylabel(a,'-log10(p)','FontSize',16,'Rotation',270);
a.FontSize = 16;
a.Label.Position(1) = 4.4;
xlabel('Region number','FontSize',16)
ylabel('Region number','FontSize',16)
set(gca,'FontSize',16)
saveas(gca,'~/Desktop/MIM_pvals.png','png')
print('~/Desktop/MIM_pvals.eps','-depsc')
close all