function fp_data2mim_sim(params)

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

%%
% define folders for saving results
DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim6/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/mim_save6/'; %folder for interim results
if ~exist(DIROUT1);mkdir(DIROUT1); end

%define which pipelines run with this configuration
if params.ip == 1 %default version 
    params.pips = [1:10 13:22];
elseif strcmp(params.ifilt(1),'c') %with champaign
    params.pips = [1:3];%FICPC1 to FIXPC3
elseif params.ip == 3 ||  params.ip == 8 
    params.pips = [1:6 9 21];
elseif params.ip == 9 || params.ip==12 
    params.pips = [1 3]; 
else % run only fixed, 90% and 99% pipelines and baseline
    params.pips = 1:9;
end


if params.ip==7 % source localization is varied
    %reload data from ip1 to keep them constant and only vary the source
    %localization
    params_save = params;
    load(sprintf('%s/mim_sensorsig/%d.mat',DIROUT1,params.iit));
    params = params_save;
    clear params_save
else
    
    if params.ip==4 || params.ip ==5 %snr/ noise mix is varied
        %reload data from ip1 to keep them constant and only vary the snr
        %or noise mix
        params_save = params;
        load(sprintf('%s/mim_sig/%d.mat',DIROUT1,params.iit));
        params= params_save;
        clear params_save
    else
        
        %% signal generation
        tic
        % getting atlas, voxel and roi indices; active voxel of each region
        % is aleady selected here
        fprintf('Getting atlas positions... \n')
        D = fp_get_Desikan(params.iReg);
        
        %signal generation
        fprintf('Signal generation... \n')
        [sig,brain_noise,sensor_noise,L,iroi_seed, iroi_tar,D, fres, n_trials,filt] = ...
            fp_generate_signal_with_timestructure(params,D,DIROUT1);
        
        if params.ip==1 %if ip1, save sig for ip4 and ip5
            dir1 =  sprintf('%s/mim_sig/',DIROUT1);
            if ~exist(dir1); mkdir(dir1); end
            outname = sprintf('%s/mim_sig/%d.mat',DIROUT1,params.iit);
            save(outname,'-v7.3')
        end
        
        t.signal = toc;
    end
    
    tic
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
    
    t.snr = toc;
    
    if params.ip==1 %if ip1, save sig for ip7
        dir1 =  sprintf('%s/mim_sensorsig/',DIROUT1);
        if ~exist(dir1); mkdir(dir1); end
        outname = sprintf('%s/mim_sensorsig/%d.mat',DIROUT1,params.iit);
        save(outname,'-v7.3')
    end
    
end

%% leadfield and inverse solution

tic

%select only voxels that belong to any roi
L_backward = L(:, D.ind_cortex, :);
ndim = size(L_backward,3); % 3 spatial dimensions 

%construct source projection filter
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
    
    [~,~,w] = awsm_champ(signal_sensor(:, :), L_perm(:, :), sigu, 200, 3, 0, 0);
    A = real(reshape(w',size(L_perm)));
    
elseif strcmp(params.ifilt,'cr') %champ with regularization 
    
    regu = fp_champ_crossval(signal_sensor,L_backward,5);
    
    sigu = regu*eye(n_sensors);
    L_perm = permute(L_backward,[1 3 2]);
    
    [~,~,w] = awsm_champ(signal_sensor(:, :), L_perm(:, :), sigu, 100, 3, 0, 0);
    A = real(reshape(w',size(L_perm)));
    
        
elseif strcmp(params.ifilt,'che') %champ heteroscedastic
    
    %estimate noise_cov
    Y = mean(signal_sensor,3);
    normalize_constant = norm(Y*Y','fro');
    [~,~,noise_cov,~,~,~,~,~,~]=nut_reg_vbfa(Y,5,150,0);
    noise_cov = noise_cov / normalize_constant;
    
    L_perm = permute(L_backward,[1 3 2]);
    
    [~,~,w_hetero,~,~,~,~]= champ_heteroscedastic_3D(signal_sensor(:, :),...
        L_perm(:, :),noise_cov,100,ndim,0,0,0,3,1);
    
    A = real(reshape(w_hetero',size(L_perm)));
    
elseif strcmp(params.ifilt,'cho') %champ homoscdastic
    
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
    
elseif strcmp(params.ifilt,'cfun') %see Hashemi 2021
    
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

t.filter = toc;


%% Loop over different pipelines

errorpipeline = [];
for ipip = params.pips 
    
    %1 to 6: FIXPC
    %7: VARPC90
    %8: VARPC99
    %9: TRUEVOX
    %10: FCMEAN
    %11: VARPC90 + correction (not in paper)
    %12: VARPC99 + correction (not in paper)
    %13 to 20: equal to 1 to 8 but with zscoring (not in paper)
    %21: MEANFC
    %22: CENTRAL
    
    fprintf(['Testing pipeline ' num2str(ipip) '\n'])
    try
        tic
        %% dimensionality reduction
        
        clear npcs  variance_explained
        signal_roi = [];
        empty_rois =[];
        active_rois = [];
        
        %loop over regions
        for aroi = 1:D.nroi
            
            clear A_ signal_source
            
            if ipip == 9 %TRUEVOX: pre-select voxels of true activity
                A_ = A(:, :,D.ind_roi_cortex{aroi},:);
                A_ = A_(:,:,D.sub_ind_roi_region{aroi},:);
                
            elseif ipip == 22 %CENTRAL: pick central voxel of each region
                A_ = A(:, :,D.ind_roi_cortex{aroi},:);
                A_ = A_(:,:,D.ctr_ind_roi_region{aroi},:);
                
            else %select all voxels of current region 
                A_ = A(:, :,D.ind_roi_cortex{aroi},:);
            end
            
            %number of voxels at the current region
            nvoxroi(aroi) = size(A_,3);
            A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
            
            %project sensor signal to voxels at the current roi (aroi)
            signal_source = A2{aroi}' * signal_sensor(:,:);
            
            if ipip > 12 && ipip < 21 %zscoring only pipelines 13 to 20 
                signal_source = zscore(signal_source')';
            end
            
            clear s_
            if strcmp(params.ifilt(1),'c') %in case of the champagne filter
                %sort out the dims and voxels with zero activity
                s_ = sum(signal_source,2)>10^-8 | sum(signal_source,2)< -(10^-8) ;
            else
                s_ = logical(ones(size(signal_source,1),1));
            end
            
            if (ipip == 9 && params.ip~=3)|| ipip == 22                
                %baseline and central voxel pipeline
                signal_roi = cat(1,signal_roi,reshape(signal_source,[],l_epoch,n_trials));
                npcs(aroi) = nvoxroi(aroi)*ndim;
 
            elseif ~(ipip == 10) && ipip < 21
                
                if strcmp(params.dimred,'p')
                    %do PCA
                    clear signal_roi_ S
                    [signal_roi_,S,~] = svd(double(signal_source(s_,:))','econ');
                    
                elseif strcmp(params.dimred,'s')
                    %do SSD (not in paper, does not work well)                   
                    [W, ~, ~, ~] = ssd(double(signal_source(s_,:))', [8 13; 2 40; 7 14], fres, [], []);
                    
                    npcs(aroi) = ipip; %fixed number of pcs
                    signal_roi_ = double(signal_source(s_,:))' * W(:, 1:npcs(aroi));
     
                else
                    error('Dimred has to be either p or s.')
                end
                
                if params.ip ~= 8 %for all pca pipelines
                    
                    % variance explained
                    vx_ = cumsum(diag(S).^2)./sum(diag(S).^2);
                    invx = 1:min(length(vx_), n_sensors);
                    varex = vx_(invx);

                    %save npcs and variance explained
                    if ismember(ipip,[7 11 19])
                        %npcs are selected in a way that 90% of the variance is preserved
                        npcs(aroi) = min(find(varex> 0.9));
                        var_explained=0.9;
                        
                    elseif ismember(ipip,[8 12 20])
                        %npcs are selected in a way that 99% of the variance is preserved
                        try %might be empty in case of the cr filter
                            npcs(aroi) = min(find(varex> 0.99));
                        catch
                            npcs(aroi) = 0;
                        end
                        var_explained=0.99;
                        
                    elseif ipip <= 6
                        %fixed number of pcs
                        npcs(aroi) = ipip;
                        try %may not be possible with the champaign filter
                            var_explained(aroi) = varex(npcs(aroi));
                        end
                        
                    elseif ipip >= 13 && ipip <= 18
                        %fixed number of pcs
                        npcs(aroi) = ipip-12;
                        var_explained(aroi) = varex(npcs(aroi));
                        
                    elseif ipip==9 && params.ip == 3
                        npcs(aroi) = nvoxroi(aroi)*ndim;
                        try 
                            var_explained(aroi) = varex(npcs(aroi));
                        end
                    end
                end
                
                if strcmp(params.ifilt(1),'c') %champaign filter
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
                
            elseif ipip == 10
                %FCMEAN
                signal_roi{aroi} = reshape(signal_source,[],l_epoch,n_trials);
                
            elseif ipip == 21
                %MEANFC
                signal_roi = cat(1,signal_roi, squeeze(mean(reshape(signal_source,ndim, nvoxroi(aroi),l_epoch,n_trials),2)));
                npcs(aroi) = ndim;
            end
        end
        
        
        
        %% calculate connectivity scores
        
        if ipip ~= 10 && ipip ~= 11 && ipip ~= 12
            
            if strcmp(params.ifilt(1),'c')
                %no calculations for empty rois
                npcs(empty_rois) = [];
            end
            
            % calculate indices
            clear inds PCA_inds
            [inds, PCA_inds] = fp_npcs2inds(npcs);
            
            % calculate FC scores 
            output = {'MIM','MIC','GC','TRGC','COH'};
            conn = data2sctrgcmim(signal_roi, fres,30, 0,0, [], inds, output,0);
            
            % extract measures out of the conn struct
            [MIM_, MIC_, GC_, DIFFGC_, iCOH_, aCOH_] = fp_unwrap_conn(conn,numel(npcs),filt,PCA_inds);
            
            if strcmp(params.ifilt(1),'c')
                % fill active rois again into the full nroixnroi structure
                
                clear mim mic diffgc
                mim = zeros(D.nroi,D.nroi);
                mim(active_rois,active_rois) = MIM_;
                MIM_ = mim;
                
                mic = zeros(D.nroi,D.nroi);
                mic(active_rois,active_rois) = MIC_;
                MIC_ = mic;
                
                gc = zeros(D.nroi,D.nroi);
                gc(active_rois,active_rois) = GC_;
                GC_ = gc;
                
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
            
        elseif ipip == 10 %FCMEAN pipeline
            
            [MIM_, MIC_] = fp_sumVox_pip(signal_roi, D.nroi, ndim, nvoxroi, fres, filt);
            
        elseif ipip == 11 || ipip == 12 %normalized VARPC pipelines (not in paper)
            
            % calculate indices
            clear inds PCA_inds
            [inds, PCA_inds] = fp_npcs2inds(npcs);
            
            %calculate normalized values of MIM and MIC
            [MIM_, MIC_, iCOH_, aCOH_] = fp_correct_mim(signal_roi,inds,fres,D,filt,PCA_inds);
            
        end
        
        t.pips(ipip) = toc;
        
        %% save stuff
        
        MIM{ipip} = MIM_;
        MIC{ipip} = MIC_;
        
        if ipip ~= 10
            
            aCOH{ipip} = aCOH_;
            iCOH{ipip} = iCOH_;
            
            if ipip ~= 11 && ipip ~= 12 
                DIFFGC{ipip} = DIFFGC_;
                GC{ipip} = GC_;
            end
            
            if ipip ~= 9 
                try
                    to_save{ipip}.npcs = npcs;
                    to_save{ipip}.varex = var_explained;
                    
                    %calculate some correlations 
                    nvoxroi_all = nvoxroi'*nvoxroi;
                    nvoxroi_all = nvoxroi_all(:);
                    corr_voxmim(ipip) = corr(nvoxroi_all,MIM_(:));
                    corr_voxmic(ipip) = corr(nvoxroi_all ,MIC_(:));
                    corr_voxicoh(ipip) = corr(nvoxroi_all,iCOH_(:));
                    corr_voxacoh(ipip) = corr(nvoxroi_all,aCOH_(:));
                    if ~strcmp(params.ifilt(1),'c')
                        corr_voxnpcs(ipip) = corr(nvoxroi', npcs');
                    end
                end
            end
            
        end
        
        
        
        %% Evaluate results with the percentile rank 
        
        [pr_mic(ipip)] = fp_pr(MIC_,iroi_seed,iroi_tar,1);        
        [pr_mim(ipip)] = fp_pr(MIM_,iroi_seed,iroi_tar,1);        
                        
        if ipip ~= 10
            [pr_aCoh(ipip)] = fp_pr(aCOH_,iroi_seed,iroi_tar,1);                  
            [pr_iCoh(ipip)] = fp_pr(iCOH_,iroi_seed,iroi_tar,1);
                   
            if ipip ~= 11 && ipip ~= 12 
                %absolute value of gc and only triu is considered. Metric neglects
                %the direction of the interaction
                [pr_abstrgc(ipip)] = fp_pr(abs(DIFFGC_),iroi_seed,iroi_tar,1);
                [pr_absgc(ipip)] = fp_pr(abs(GC_),iroi_seed,iroi_tar,1);
                
                %only positive part of gc is submitted and the whole matrix is
                %considered. Metric that is strongly influenced by the direction of
                %the effect
                clear pos_difftrgc
                pos_difftrgc = DIFFGC_;
                pos_difftrgc(pos_difftrgc< 0) = 0;
                [pr_postrgc(ipip)]...
                    = fp_pr(pos_difftrgc,iroi_seed,iroi_tar,0);
                
                clear pos_diffgc
                pos_diffgc = GC_;
                pos_diffgc(pos_diffgc< 0) = 0;
                [pr_posgc(ipip)]...
                    = fp_pr(pos_diffgc,iroi_seed,iroi_tar,0);
                
            end           
        end
                     
        clear MIM_ MIC_ DIFFGC_ GC_ aCOH_ iCOH_
        
    catch
        errorpipeline = [errorpipeline ipip]
    end
    
end %pips

%% Saving
fprintf('Saving... \n')
%save all
outname = sprintf('%smim_%s.mat',DIROUT,params.logname);
save(outname,'-v7.3')

%save only evaluation parameters
outname1 = sprintf('%spr_%s.mat',DIROUT,params.logname);
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

