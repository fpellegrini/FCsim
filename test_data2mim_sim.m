%% parameter configuration

DIROUT1=[];
params.iReg=1;
params.iInt = 2;
params.ilag = 2;
params.isnr = 0.7;
params.iss = 0.5;
params.ip=2;
params.ihemi = 0;
zs = 0;
params.ifilt = 'l';
t=[];

%%
for iit = 1:5
    %% signal generation
    
    D = fp_get_Desikan(params.iReg);
    
    [sig,brain_noise,sensor_noise,L,iroi_seed, iroi_tar,D, fres, n_trials,filt] = ...
        fp_generate_signal_with_timestructure(params,D,DIROUT1);
    
    
    %combine noise sources
    noise = params.iss*brain_noise + (1-params.iss)*sensor_noise;
    noise_f = (filtfilt(filt.bband, filt.aband, noise'))';
    noise = noise ./ norm(noise_f, 'fro');
    %combine signal and noise
    signal_sensor1 = params.isnr*sig + (1-params.isnr)*noise;
    signal_sensor_f = (filtfilt(filt.bband, filt.aband, signal_sensor1'))';
    signal_sensor1 = signal_sensor1 ./ norm(signal_sensor_f, 'fro');
    
    signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
    signal_sensor = signal_sensor / norm(signal_sensor, 'fro');
    
    %reshape
    signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
    [n_sensors, l_epoch, n_trials] = size(signal_sensor);
    
    
    tic
    
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
        
    elseif strcmp(params.ifilt,'c') %champ
        sigma2_total = mean(diag(cov(signal_sensor(:, :)')));
        regu = sigma2_total*0.2;
        sigu = regu*eye(n_sensors);
        L_perm = permute(L_backward,[1 3 2]);
        
        [~,~,w] = awsm_champ(signal_sensor(:, :), L_perm(:, :), sigu, 200, 3, 2, 0);
        A = real(reshape(w',size(L_perm)));
    end
    
    t.filter = toc;
    
    
    %% Loop over different pipelines
    
    
    %%%%%%%to do: in most of the cases only pips 1 to 9%%%%%%%%%%%%
    for ipip = 3
        
        %1 to 6: fixed
        %7: 90%
        %8: 99%
        %9: baseline
        %10: sumVox
        %11: 90% corrected
        %12: 99% corrected
        %13 to 20: equal to 1 to 8 but with zscoring
        
        fprintf(['Testing pipeline ' num2str(ipip) '\n'])
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
            
            if ipip == 9 %baseline: pre-select voxels of true activity
                A_ = A_(:,:,D.sub_ind_roi_region{aroi});
                
            elseif ipip == 22 %pick central voxel of each region
                A_ = A_(:,:,D.ctr_ind_roi_region{aroi});
            end
            
            %number of voxels at the current roi
            nvoxroi(aroi) = size(A_,3);
            
            A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
            
            %project sensor signal to voxels at the current roi (aroi)
            signal_source = A2{aroi}' * signal_sensor(:,:);
            
            if ipip > 12 %zscoring pipelines
                signal_source = zscore(signal_source);
            end
            
            if strcmp(params.ifilt,'c')
                %sort out the dims and voxels with zero activity
                s_ = sum(signal_source,2)>10^-8 | sum(signal_source,2)< -(10^-8) ;
            else
                s_ = logical(ones(size(signal_source,1),1));
            end
            
            if ~(ipip == 9) && ~(ipip == 10) && ipip < 21
                
                %do PCA
                clear signal_roi_ S
                [signal_roi_,S,~] = svd(double(signal_source(s_,:))','econ');
                
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
                    npcs(aroi) = min(find(varex> 0.99));
                    var_explained=0.99;
                elseif ipip <= 6
                    %fixed number of pcs
                    npcs(aroi) = ipip;
                    try
                        var_explained(aroi) = varex(npcs(aroi));
                    end
                elseif ipip >= 13 && ipip <= 18
                    %fixed number of pcs
                    npcs(aroi) = ipip-12;
                    var_explained(aroi) = varex(npcs(aroi));
                end
                
                if strcmp(params.ifilt,'c')
                    %sort out voxels with zero activity
                    try
                        %bring signal_roi to the shape of npcs x l_epoch x n_trials
                        signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
                        active_rois = [active_rois aroi];
                    catch
                        empty_rois = [empty_rois aroi];
                    end
                    
                else
                    %bring signal_roi to the shape of npcs x l_epoch x n_trials
                    signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
                end
                
            elseif ipip == 9 || ipip == 22
                %baseline
                signal_roi = cat(1,signal_roi,reshape(signal_source,[],l_epoch,n_trials));
                npcs(aroi) = 3;
                
            elseif ipip == 10
                %sumVox
                signal_roi{aroi} = reshape(signal_source,[],l_epoch,n_trials);
                
            elseif ipip == 21
                
                signal_roi = cat(1,signal_roi, squeeze(mean(reshape(signal_source,ndim, nvoxroi(aroi),l_epoch,n_trials),2)));
                npcs(aroi) = 3;
            end
        end
        
        
        
        %% calculate connectivity scores
        
        if ipip ~= 10 && ipip ~= 11 && ipip ~= 12
            
            if strcmp(params.ifilt,'c')
                npcs(empty_rois) = [];
            end
            
            % calculate indices
            clear inds PCA_inds
            [inds, PCA_inds] = fp_npcs2inds(npcs);
            
            % calculate MIM, MIC, TRGC and COH
            output = {'MIM','MIC','COH'};
            conn = data2sctrgcmim(signal_roi, fres, 20, 0,0, [], inds, output);
            
            % extract measures out of the conn struct
            [MIM_, MIC_, DIFFGC_, iCOH_, aCOH_] = fp_unwrap_conn(conn,numel(npcs),filt,PCA_inds);
            
            
            if strcmp(params.ifilt,'c') && ~(isempty(empty_rois))
                clear mim mic diffgc
                mim = zeros(D.nroi,D.nroi);
                mim(active_rois,active_rois) = MIM_;
                MIM_ = mim;
                
                mic = zeros(D.nroi,D.nroi);
                mic(active_rois,active_rois) = MIC_;
                MIC_ = mic;
                
%                 diffgc = zeros(D.nroi,D.nroi);
%                 diffgc(active_rois,active_rois) = DIFFGC_;
%                 DIFFGC_ = diffgc;
                
                icoh = zeros(D.nroi,D.nroi);
                icoh(active_rois,active_rois) = iCOH_;
                iCOH_ = icoh;
                
                acoh = zeros(D.nroi,D.nroi);
                acoh(active_rois,active_rois) = aCOH_;
                aCOH_ = acoh;                
            end
            
        elseif ipip == 10 %sumVox pipeline
            
            [MIM_, MIC_] = fp_sumVox_pip(signal_roi, D.nroi, ndim, nvoxroi, fres, filt);
            
        elseif ipip == 11 || ipip == 12 %normalized 90% and 99% pipeline
            
            % calculate indices
            clear inds PCA_inds
            [inds, PCA_inds] = fp_npcs2inds(npcs);
            
            %calculate normalized values of MIM and MIC
            [MIM_, MIC_, iCOH_, aCOH_] = fp_correct_mim(signal_roi,inds,fres,D,filt,PCA_inds);
            
        end
        
        
        t.pips(ipip) = toc;
        
        %% save a few more things
        
        
        [mrr_mic(ipip,iit), pr_mic(ipip),hk_mic(ipip)] = fp_mrr_hk(MIC_,iroi_seed,iroi_tar,1);
        [mrr_mim(ipip,iit), pr_mim(ipip),hk_mim(ipip)] = fp_mrr_hk(MIM_,iroi_seed,iroi_tar,1);
        
        if ipip ~= 10
            [mrr_aCoh(ipip,iit), pr_aCoh(ipip),hk_aCoh(ipip)] = fp_mrr_hk(aCOH_,iroi_seed,iroi_tar,1);
            [mrr_iCoh(ipip,iit), pr_iCoh(ipip),hk_iCoh(ipip)] = fp_mrr_hk(iCOH_,iroi_seed,iroi_tar,1);
            
%             %absolute value of gc and only triu is considered. Metric neglects
%             %the direction of the interaction
%             [mrr_absgc(ipip,iit), pr_absgc(ipip),hk_absgc(ipip)] = fp_mrr_hk(abs(DIFFGC_),iroi_seed,iroi_tar,1);
%             
%             %only positive part of gc is submitted and the whole matrix is
%             %considered. Metric that is strongly influenced by the direction of
%             %the effect
%             pos_diffgc = DIFFGC_;
%             pos_diffgc(pos_diffgc < 0) = 0;
%             [mrr_posgc(ipip,iit), pr_posgc(ipip),hk_posgc(ipip)] = fp_mrr_hk(pos_diffgc,iroi_seed,iroi_tar,0);

        end
        
        
        clear MIM_ MIC_ DIFFGC_ aCOH_ iCOH_
        
        
    end %pips
    
end

%% plot

for ipip = 3
    
    figure
    subplot(2,2,1)
    fp_raincloud_plot(mrr_mic(ipip,:), [0.7 0.7 0.8], 1,0.2, 'ks');
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    set(gca, 'XLim', [0 1]);
    title('MIC')
    
    subplot(2,2,2)
    fp_raincloud_plot(mrr_mim(ipip,:), [0.7 0.7 0.8], 1,0.2, 'ks');
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    set(gca, 'XLim', [0 1]);
    title('MIM')
    
    subplot(2,2,3)
    fp_raincloud_plot(mrr_iCoh(ipip,:), [0.7 0.7 0.8], 1,0.2, 'ks');
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    set(gca, 'XLim', [0 1]);
    title('iCOH')
    
    subplot(2,2,4)
    fp_raincloud_plot(mrr_aCoh(ipip,:), [0.7 0.7 0.8], 1,0.2, 'ks');
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    set(gca, 'XLim', [0 1]);
    title('aCOH')
    
end



