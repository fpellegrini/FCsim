function fp_lag_sim(params)

% define folders for saving results
DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim5_lag/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/mim_save5_lag/';
if ~exist(DIROUT1);mkdir(DIROUT1); end

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
        
    elseif strcmp(params.ifilt,'l') %lcmv (default)
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
        
        [~,~,w] = awsm_champ(signal_sensor(:, :), L_perm(:, :), sigu, 100, 3, 2, 0);
        A = real(reshape(w',size(L_perm)));
        
        
    elseif strcmp(params.ifilt,'che') %champ hetero
        
        %estimate noise_cov
        Y = mean(signal_sensor,3);
        normalize_constant = norm(Y*Y','fro');
        [~,~,noise_cov,~,~,~,~,~,~]=nut_reg_vbfa(Y,5,150,0);
        noise_cov = noise_cov / normalize_constant;
        
        L_perm = permute(L_backward,[1 3 2]);
        
        [~,~,w_hetero,~,~,~,~]= champ_heteroscedastic_3D(signal_sensor(:, :),...
            L_perm(:, :),noise_cov,100,ndim,0,0,0,3,1);
        
        A = real(reshape(w_hetero',size(L_perm)));
        
    elseif strcmp(params.ifilt,'cho') %chanp homo
        
        %estimate noise_cov
        Y = mean(signal_sensor,3);
        normalize_constant = norm(Y*Y','fro');
        [~,~,noise_cov,~,~,~,~,~,~]=nut_reg_vbfa(Y,5,150,0);
        noise_cov = eye(n_sensors) * mean(diag(noise_cov));
        noise_cov = noise_cov / normalize_constant;
        
        L_perm = permute(L_backward,[1 3 2]);
        L_perm = L_perm / normalize_constant;
        signal_sensor_n = signal_sensor./normalize_constant;
        
        [~,~,~,~,W_Champ_Homoscedastic]= Champagne_Homoscedastic(L_perm(:,:),signal_sensor_n(:,:),...
            'noise_cov',noise_cov,'noise_update_mode',1,'max_num_iterations',100,...
            'threshold_error',eps,'print_results',0,'print_figures',0);
        
        A = real(reshape(W_Champ_Homoscedastic',size(L_perm)));
        
    elseif strcmp(params.ifilt,'cfun')
        
        %estimate noise_cov
        Y = mean(signal_sensor,3);
        normalize_constant = norm(Y*Y','fro');
        [~,~,noise_cov,~,~,~,~,~,~]=nut_reg_vbfa(Y,5,150,0);
        noise_cov = noise_cov / normalize_constant;
        
        L_perm = permute(L_backward,[1 3 2]);
        
        [~,~,~,~,W_FUN]= FUN_learning_func_3D(L_perm(:,:),signal_sensor(:,:),...
            'noise_cov',noise_cov,'noise_update_mode','Diagonal',...
            'source_update_mode','Diagonal','max_num_iterations',100,...
            'threshold_error',eps,'print_results',0,'print_figures',0);
        
        
        A = real(reshape(W_FUN',size(L_perm)));
        
    end
    
    %% Pipeline 1
    
    errorpipeline = [];
    for ipip = 3 %fixed 3 pc
        
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
               
                %zscoring
                signal_source = zscore(signal_source')';
                               
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
            
            output = {'MIM','MIC','GC','TRGC','COH'};
            conn = data2sctrgcmim(signal_roi, fres,30, 0,0, [], inds, output,0);
            
            % extract measures out of the conn struct
            [MIM_, MIC_, GC_, DIFFGC_, iCOH_, aCOH_] = fp_unwrap_conn(conn,numel(npcs),filt,PCA_inds);
            
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
            GC{ipip} = GC_;
            
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
        [pr_mic(ipip)] = fp_mrr_hk_short(MIC_,iroi_seed,iroi_tar,1);        
        [pr_mim(ipip)] = fp_mrr_hk_short(MIM_,iroi_seed,iroi_tar,1);        
                
        
        if ipip ~= 10
            [pr_aCoh(ipip)] = fp_mrr_hk_short(aCOH_,iroi_seed,iroi_tar,1);                  
            [pr_iCoh(ipip)] = fp_mrr_hk_short(iCOH_,iroi_seed,iroi_tar,1);
                   
            if ipip ~= 11 && ipip ~= 12 
                %absolute value of gc and only triu is considered. Metric neglects
                %the direction of the interaction
                [pr_abstrgc(ipip)] = fp_mrr_hk_short(abs(DIFFGC_),iroi_seed,iroi_tar,1);
                [pr_absgc(ipip)] = fp_mrr_hk_short(abs(GC_),iroi_seed,iroi_tar,1);
                
                %only positive part of gc is submitted and the whole matrix is
                %considered. Metric that is strongly influenced by the direction of
                %the effect
                clear pos_difftrgc
                pos_difftrgc = DIFFGC_;
                pos_difftrgc(pos_difftrgc< 0) = 0;
                [pr_postrgc(ipip)]...
                    = fp_mrr_hk_short(pos_difftrgc,iroi_seed,iroi_tar,0);
                
                clear pos_diffgc
                pos_diffgc = GC_;
                pos_diffgc(pos_diffgc< 0) = 0;
                [pr_posgc(ipip)]...
                    = fp_mrr_hk_short(pos_diffgc,iroi_seed,iroi_tar,0);
                
            end
            
        end
              
        
        clear MIM_ MIC_ DIFFGC_ GC_ aCOH_ iCOH_
        
    end %pips
    
    %% Saving
    fprintf('Saving... \n')
    %save all
    outname = sprintf('%smim_%s_lag%s.mat',DIROUT,params.logname,num2str(ilag));
    save(outname,'-v7.3')
    
    %save only evaluation parameters
    outname1 = sprintf('%spr_%s_lag%s.mat',DIROUT,params.logname,num2str(ilag));
    save(outname1,...
    'pr_mic',...
    'pr_mim',...
    'pr_aCoh',...
    'pr_iCoh',...
    'pr_absgc',...
    'pr_posgc',...
    'pr_abstrgc',...
    'pr_postrgc',...
    '-v7.3')
    
end

