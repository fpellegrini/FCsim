function fp_simulation_stats

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

rng('default')
rng(5)

params.iReg=1; %number of interacting voxels in interacting regions 
params.iInt = 1; %number of interacting regions
params.ilag = 2; %lag size
params.isnr = 0.6; %SNR
params.iss = 0.5; %noise mix
t=[]; %run time
params.ifilt='l';
params.dimred='p';
params.iit = 1;
params.ip = NaN;
ipip=3;

                            
%% signal generation

% getting atlas, voxel and roi indices; active voxel of each region
% is aleady selected here
fprintf('Getting atlas positions... \n')
D = fp_get_Desikan(params.iReg);

%signal generation
fprintf('Signal generation... \n')
[sig,brain_noise,sensor_noise,L,iroi_seed, iroi_tar,D, fres, n_trials,filt,~] = ...
    fp_u(params,D,[]);

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


%% leadfield and inverse solution

%select only voxels that belong to any roi
L_backward = L(:, D.ind_cortex, :);
ndim = size(L_backward,3); % 3 spatial dimensions

cCS = cov(signal_sensor(:,:)');
reg = 0.05*trace(cCS)/length(cCS);
Cr = cCS + reg*eye(size(cCS,1));

[~, A] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
A = permute(A,[1, 3, 2]);


%% PCA

clear npcs  variance_explained
signal_roi = [];
empty_rois =[];
active_rois = [];

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
    try %may not be possible with the champaign filter
        var_explained(aroi) = varex(npcs(aroi));
    end
    
    %bring signal_roi to the shape of npcs x l_epoch x n_trials
    signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
    
end

%% calculate true MIM

[nchan, ~, nepo] = size(signal_roi);
maxfreq = fres+1;

% calculate indices
clear inds PCA_inds
[inds, PCA_inds] = fp_npcs2inds(npcs);
ninds = length(inds);

%true MIM 
output = {'MIM'};
conn = data2sctrgcmim(signal_roi, fres, 30, 0,0, [], inds, output);
[MIM_t, ~, ~, ~, ~,~] = fp_unwrap_conn(conn,D.nroi,filt,PCA_inds);

outname1 = ['~/Desktop/signal.mat'];
save(outname1,'signal_roi','iroi_seed','iroi_tar','npcs','MIM_t','D','filt','-v7.3')

