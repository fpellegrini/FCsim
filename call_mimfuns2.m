
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

fres = 40;
n_trials = 200;

%% signal generation

% ROI labels
% In D.sub_ind_roi, there are the randomly
% selected voxels of each region
fprintf('Getting atlas positions... \n')
D = fp_get_Desikan(params.iReg);

%signal generation
fprintf('Signal generation... \n')
[sig,brain_noise,sensor_noise,gt,L,iroi_seed, iroi_tar,D] = fp_generate_mim_signal(params, ...
    fres,n_trials, D,[]);

%combine noise sources
for itrial = 1:n_trials
    noise = params.iss*brain_noise{itrial} + (1-params.iss)*sensor_noise{itrial};
    noise = noise ./ norm(noise, 'fro');
    %combine signal and noise
    signal_sensor1 = params.isnr*sig{itrial} + (1-params.isnr)*noise;
    signal_sensor(:,:,itrial) = signal_sensor1 ./ norm(signal_sensor1, 'fro');
end


%% get CS and filter A
%parameters
id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
id_meg_chan = 1:size(signal_sensor,1);
nmeg = numel(id_meg_chan);

%cross spectrum
fprintf('Calculating cross spectrum... \n')
CS = fp_tsdata_to_cpsd(signal_sensor,fres,'WELCH',...
    [id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
CS(:,:,1)=[];
nfreq = size(CS,3);

% leadfield
L3 = L(:, D.ind_cortex, :);
L_backward = L3;
ni = size(L_backward,3);

%construct source filter
if strcmp(params.ifilt,'e')
    A = squeeze(mkfilt_eloreta_v2(L_backward));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(params.ifilt,'d')
    
    A=zeros(nmeg,ni,D.nvox,nfreq);
    
    for ifrq = 1:nfreq
        cCS = CS(:,:,ifrq);
        lambda = mean(diag(real(cCS)))/100;
        
        CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
        
        for is=1:D.nvox %iterate across nodes
            Lloc=squeeze(L_backward(:,is,:));
            A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
        end
    end
    fqA = 1:nfreq; %This filter is frequency specific.
    nfqA = nfreq;
    
    
elseif strcmp(params.ifilt,'l')
    cCS = sum(CS,3);
    reg = 0.05*trace(cCS)/length(cCS);
    Cr = cCS + reg*eye(size(cCS,1));
    
    [~, A] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
end

%% calculate filter and npcs

mode1 = 'max'; % I took out the option to call 'all' at once here for a better overview 

for aroi = 1:D.nroi
    
    %filter at current roi
    clear A_ CSv
    A_ = A(:, :,D.ind_roi_cortex{aroi},:);
    nvoxroi = size(A_,3); %voxels in the current roi
    A2{aroi} = reshape(A_, [nmeg, ni*nvoxroi, nfqA]);
    
    
    if ~strcmp(mode1,'case2')&& ~strcmp(mode1,'baseline')&& ~strcmp(mode1,'bandc')
        
        %project CS to voxel space
        for ifq = 1: nfreq
            CSv(:,:,ifq) = squeeze(A2{aroi}(:,:,fqA(ifq)))' * CS(:,:,ifq)...
                * squeeze(A2{aroi}(:,:,fqA(ifq)));
        end
        
        %zscoring
        clear CSz
        ZS{aroi} = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
        
        %apply ZS
        for ifreq = 1:nfreq
            CSz(:,:, ifreq) = ZS{aroi}'*squeeze(CSv(:, :,ifreq))*ZS{aroi};
        end
        
        %PCA
        clear CSs v v5 in V_ D_
        CSs = squeeze(sum(real(CSz),3)); %covariance
        [V_, D_] = eig(CSs);
        [D_, in] = sort(real(diag(D_)), 'descend');
        V{aroi} = V_(:,in);
        
        %npcs
        if isnumeric(mode1)
            %fixed number of pcs for every roi
            npcs(aroi) = mode1;
            
        elseif strcmp(mode1,'max')
            %pipeline 6)
            
            vx_ = cumsum(D_)./sum(D_);
            npcs(aroi) = min(find(vx_>0.99999));
            
        elseif strcmp(mode1,'percent')
            %pipeline 7)
            
            % variance explained
            vx_ = cumsum(D_)./sum(D_);
            npcs(aroi) = min(find(vx_>0.9));            
        end
    else
        % if mode1 is case2 and/or baseline 
        ZS = [];
        V=[];
        npcs=[];
    end
end

%% calculate MIM/MIC

[mic,mim,to_save] = fp_compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA,nfqA,params.ihemi);




