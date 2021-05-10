function fp_lag_sim(params)

% define folders for saving results
DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim4_lag/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/mim_save4_lag/';
if ~exist(DIROUT1);mkdir(DIROUT1); end

%define which pipelines run with this configuration
params.pips = 1;

%% signal generation
% getting atlas, voxel and roi indices; active voxel of each region
% is aleady selected here
fprintf('Getting atlas positions... \n')
D = fp_get_Desikan(params.iReg);

%% vary lags

for ilag = 1:5 %in ms: 2,4,6,8,10
    
    %signal generation
    fprintf('Signal generation... \n')
    [sig,brain_noise,sensor_noise,L,iroi_seed, iroi_tar,D, fres, n_trials,filt] = ...
        fp_generate_signal_lagsim(params,D,DIROUT1,ilag);
    
    
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
    ndim = size(L_backward,3);
    
    %construct source filter
    if strcmp(params.ifilt,'e') %eloreta
        reg_param = fp_eloreta_crossval(signal_sensor,L_backward,5);
        A = squeeze(mkfilt_eloreta_v2(L_backward,reg_param));
        A = permute(A,[1, 3, 2]);
        
    elseif strcmp(params.ifilt,'l') %lcmv
        cCS = cov(signal_sensor(:,:)');
        reg = 0.05*trace(cCS)/length(cCS);
        Cr = cCS + reg*eye(size(cCS,1));
        
        [~, A] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
        A = permute(A,[1, 3, 2]);
        
    elseif strcmp(params.ifilt,'c') %champ with non-realistically good conditions
        sigma2_total = mean(diag(cov(signal_sensor(:, :)')));
        regu = sigma2_total* (1-params.isnr);
        
        sigu = regu*eye(n_sensors);
        L_perm = permute(L_backward,[1 3 2]);
        
        [~,~,w] = awsm_champ(signal_sensor(:, :), L_perm(:, :), sigu, 200, 3, 2, 0);
        A = real(reshape(w',size(L_perm)));
        
    elseif strcmp(params.ifilt,'cr') %champ with regulaization
        
        regu = fp_champ_crossval(signal_sensor,L_backward,5);
        
        sigu = regu*eye(n_sensors);
        L_perm = permute(L_backward,[1 3 2]);
        
        [~,~,w] = awsm_champ(signal_sensor(:, :), L_perm(:, :), sigu, 200, 3, 2, 0);
        A = real(reshape(w',size(L_perm)));
        
    end
    
    %% Pipeline 1
    
    errorpipeline = [];
    for ipip = 1 %fixed 1 pc
        
        fprintf(['Testing pipeline ' num2str(ipip) '\n'])
        try
            tic
            %% PCA
            
            clear npcs
            signal_roi = [];
            empty_rois =[];
            active_rois = [];
            
            %loop over regions
            for aroi = 1:D.nroi
                
                clear A_ signal_source
                
                %A_ is the lcmv filter at aroi
                A_ = A(:, :,D.ind_roi_cortex{aroi},:);
                
                %number of voxels at the current roi
                nvoxroi(aroi) = size(A_,3);
                A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
                
                %project sensor signal to voxels at the current roi (aroi)
                signal_source = A2{aroi}' * signal_sensor(:,:);
                
                clear s_
                if strcmp(params.ifilt,'c') || strcmp(params.ifilt,'cr') %in case of the champaign filter
                    %sort out the dims and voxels with zero activity
                    s_ = sum(signal_source,2)>10^-8 | sum(signal_source,2)< -(10^-8) ;
                else
                    s_ = logical(ones(size(signal_source,1),1));
                end
                
                
                %do PCA
                clear signal_roi_ S
                [signal_roi_,S,~] = svd(double(signal_source(s_,:))','econ');
                
                % variance explained
                vx_ = cumsum(diag(S).^2)./sum(diag(S).^2);
                invx = 1:min(length(vx_), n_sensors);
                varex = vx_(invx);
                
                %fixed number of pcs
                npcs(aroi) = ipip;
                try %may not be possible with the champaign filter
                    var_explained(aroi) = varex(npcs(aroi));
                end
                
                if strcmp(params.ifilt,'c') || strcmp(params.ifilt,'cr') %champaign filter
                    %sort out rois with zero activity
                    if npcs(aroi) == 0
                        empty_rois = [empty_rois aroi];
                    else
                        try
                            %bring signal_roi to the shape of npcs x l_epoch x n_trials
                            signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
                            active_rois = [active_rois aroi];
                        catch
                            empty_rois = [empty_rois aroi];
                        end
                    end
                    
                else
                    %bring signal_roi to the shape of npcs x l_epoch x n_trials
                    signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
                end
                
            end
            
            
            %% calculate connectivity scores
            
            if strcmp(params.ifilt,'c') || strcmp(params.ifilt,'cr')
                %no calculations for empty rois
                npcs(empty_rois) = [];
            end
            
            % calculate indices
            clear inds PCA_inds
            [inds, PCA_inds] = fp_npcs2inds(npcs);
            
            output = {'MIM','MIC','TRGC','COH'};
            conn = data2sctrgcmim(signal_roi, fres,30, 0,0, [], inds, output,0);
            
            % extract measures out of the conn struct
            [MIM_, MIC_, DIFFGC_, iCOH_, aCOH_] = fp_unwrap_conn(conn,numel(npcs),filt,PCA_inds);
            
            if strcmp(params.ifilt,'c') || strcmp(params.ifilt,'cr')
                % fill active rois again into the full nroixnroi structure
                
                clear mim mic diffgc
                mim = zeros(D.nroi,D.nroi);
                mim(active_rois,active_rois) = MIM_;
                MIM_ = mim;
                
                mic = zeros(D.nroi,D.nroi);
                mic(active_rois,active_rois) = MIC_;
                MIC_ = mic;
                
                diffgc = zeros(D.nroi,D.nroi);
                diffgc(active_rois,active_rois) = DIFFGC_;
                DIFFGC_ = diffgc;
                
                icoh = zeros(D.nroi,D.nroi);
                icoh(active_rois,active_rois) = iCOH_;
                iCOH_ = icoh;
                
                acoh = zeros(D.nroi,D.nroi);
                acoh(active_rois,active_rois) = aCOH_;
                aCOH_ = acoh;
                
            end
            
            
            
            %% save stuff
            
            MIM{ipip} = MIM_;
            MIC{ipip} = MIC_;
            
            
            aCOH{ipip} = aCOH_;
            iCOH{ipip} = iCOH_;
            
            DIFFGC{ipip} = DIFFGC_;
            
            try
                to_save{ipip}.npcs = npcs;
                to_save{ipip}.varex = var_explained;
                
                nvoxroi_all = nvoxroi'*nvoxroi;
                nvoxroi_all = nvoxroi_all(:);
                corr_voxmim(ipip) = corr(nvoxroi_all,MIM_(:));
                corr_voxmic(ipip) = corr(nvoxroi_all ,MIC_(:));
                corr_voxicoh(ipip) = corr(nvoxroi_all,iCOH_(:));
                corr_voxacoh(ipip) = corr(nvoxroi_all,aCOH_(:));
                if ~strcmp(params.ifilt,'c') || ~strcmp(params.ifilt,'cr')
                    corr_voxnpcs(ipip) = corr(nvoxroi', npcs');
                end
            end
        end
        
        
        
        
        %% Evaluate
        
        [mrr_mic(ipip), pr_mic(ipip),hk_mic(ipip),em1_mic(ipip),em2_mic(ipip),em3_mic(ipip)] ...
            = fp_mrr_hk(MIC_,iroi_seed,iroi_tar,1);
        [mrr_mim(ipip), pr_mim(ipip),hk_mim(ipip),em1_mim(ipip),em2_mim(ipip),em3_mim(ipip)] ...
            = fp_mrr_hk(MIM_,iroi_seed,iroi_tar,1);
        
        [mrr_aCoh(ipip), pr_aCoh(ipip),hk_aCoh(ipip),em1_aCoh(ipip),em2_aCoh(ipip),em3_aCoh(ipip)] ...
            = fp_mrr_hk(aCOH_,iroi_seed,iroi_tar,1);
        [mrr_iCoh(ipip), pr_iCoh(ipip),hk_iCoh(ipip),em1_iCoh(ipip),em2_iCoh(ipip),em3_iCoh(ipip)] ...
            = fp_mrr_hk(iCOH_,iroi_seed,iroi_tar,1);
        
        %absolute value of gc and only triu is considered. Metric neglects
        %the direction of the interaction
        [mrr_absgc(ipip), pr_absgc(ipip),hk_absgc(ipip),...
            em1_absgc(ipip),em2_absgc(ipip),em3_absgc(ipip)] ...
            = fp_mrr_hk(abs(DIFFGC_),iroi_seed,iroi_tar,1);
        
        %only positive part of gc is submitted and the whole matrix is
        %considered. Metric that is strongly influenced by the direction of
        %the effect
        clear pos_diffgc
        pos_diffgc = DIFFGC_;
        pos_diffgc(pos_diffgc< 0) = 0;
        [mrr_posgc(ipip), pr_posgc(ipip),hk_posgc(ipip),...
            em1_posgc(ipip),em2_posgc(ipip),em3_posgc(ipip)] ...
            = fp_mrr_hk(pos_diffgc,iroi_seed,iroi_tar,0);
        
        %wrong directions
        clear pos_diffgc_w
        pos_diffgc_w = -DIFFGC_;
        pos_diffgc_w(pos_diffgc_w < 0) = 0;
        [mrr_posgc_w(ipip), pr_posgc_w(ipip),hk_posgc_w(ipip),...
            em1_posgc_w(ipip),em2_posgc_w(ipip),em3_posgc_w(ipip)]...
            = fp_mrr_hk(pos_diffgc_w,iroi_seed,iroi_tar,0);
        
        
        
        clear MIM_ MIC_ DIFFGC_ aCOH_ iCOH_
        
        
    end %pips
    
    %% Saving
    fprintf('Saving... \n')
    %save all
    outname = sprintf('%smim_%s_lag%s.mat',DIROUT,params.logname,num2str(ilag));
    save(outname,'-v7.3')
    
    %save only evaluation parameters
    outname1 = sprintf('%smrr_%s_lag%s.mat',DIROUT,params.logname,num2str(ilag));
    save(outname1,...
        'mrr_mic','pr_mic','hk_mic','em1_mic','em2_mic','em3_mic',...
        'mrr_mim','pr_mim','hk_mim','em1_mim','em2_mim','em3_mim',...
        'mrr_aCoh','pr_aCoh','hk_aCoh','em1_aCoh','em2_aCoh','em3_aCoh',...
        'mrr_iCoh','pr_iCoh','hk_iCoh','em1_iCoh','em2_iCoh','em3_iCoh',...
        'mrr_absgc','pr_absgc','hk_absgc','em1_absgc','em2_absgc','em3_absgc',...
        'mrr_posgc','pr_posgc','hk_posgc','em1_posgc','em2_posgc','em3_posgc',...
        'mrr_posgc_w','pr_posgc_w','hk_posgc_w','em1_posgc_w','em2_posgc_w','em3_posgc_w',...
        '-v7.3')
    
end

