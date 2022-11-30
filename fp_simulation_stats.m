fp_addpath

DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim5_stats/';
if ~exist(DIROUT);mkdir(DIROUT); end

params.pips = 3;
ipip = 3;
[iInt,iReg,isnr,iss,ilag,ifilt, dimred] = fp_get_params(3);

params.iInt = iInt; %number of interactions
params.iReg = iReg; %number of active voxels per region
params.isnr = isnr; %signal-to-noise ratio
params.iss = iss; %noise mix
params.ilag = ilag; %magitude of time delay (ilag=2 -> time delay between 50 and 200 Hz)
params.ifilt = ifilt{1}; %source projection filter type
params.dimred = dimred; %type of dimensionality reduction (pca vs ssd, but ssd not properly investigated here)
params.iit = 1;
params.ip = NaN;

logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_iter%d'...
    ,iInt,iReg,isnr*10,iss*10, ilag,ifilt{1},dimred,1);
params.logname = logname;
                            
%% signal generation

% getting atlas, voxel and roi indices; active voxel of each region
% is aleady selected here
fprintf('Getting atlas positions... \n')
D = fp_get_Desikan(params.iReg);

%signal generation
fprintf('Signal generation... \n')
[sig,brain_noise,sensor_noise,L,iroi_seed, iroi_tar,D, fres, n_trials,filt] = ...
    fp_generate_signal_with_timestructure(params,D,[]);

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

%% calculate true MIM and do statistics 
tic
nit = 1000;
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
   

for iit = 1:nit %one iteration takes ~90 sec on my local laptop

    %shuffle trials
    shuf_inds = randperm(nepo);   
    
    clear MIM2 CS COH2
    CS = fp_tsdata_to_cpsd(signal_roi, fres, 'WELCH', 1:nchan, 1:nchan,1:nepo,shuf_inds);
    
    for ifreq = 1:maxfreq
        clear pow
        pow = real(diag(CS(:,:,ifreq)));
        COH2(:,:,ifreq) = CS(:,:,ifreq)./ sqrt(pow*pow');
    end
    
    % loop over sender/receiver combinations to compute time-reversed GC
    for iind = 1:ninds
        if ~isequal(inds{iind}{1}, inds{iind}{2})
            %ind configuration
            subset = [inds{iind}{1} inds{iind}{2}];
            subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
            
            %MIC and MIM
            [~ , MIM2(:, iind)] =  roi_mim2(COH2(subset, subset, :), subinds{1}, subinds{2});
        end
    end        
      
    % extract measures out of the conn struct
    clear conn
    conn.MIM = MIM2;
    conn.inds = inds;  
    [MIM_s(:,:,iit), ~, ~, ~, ~, ~] = fp_unwrap_conn(conn,D.nroi,filt,PCA_inds);

end

for iroi = 1:D.nroi 
    for jroi = 1:D.nroi
        MIM_p(iroi,jroi) = sum(squeeze(MIM_s(iroi,jroi,:))>MIM_t(iroi,jroi))/nit;
    end
end
t = toc; 
%% save 
outname1 = [DIROUT 'result.mat'];
save(outname1,'MIM_t','MIM_s','MIM_p','signal_roi','t','-v7.3')
