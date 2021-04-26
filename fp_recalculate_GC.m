function fp_recalculate_GC

% define folders for saving results
DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim3/';
if ~exist(DIROUT);mkdir(DIROUT); end

for iname = 1
    
    iname
    
    clearvars -except iname DIRIN
    
    %default paramenters
    nit = 100;
    iInt = 2;
    iReg=1;
    isnr=0.7;
    iss = 0.5;
    ilag=2;
    ifilt='l';
    
    if iname ==  2
        iInt = 1;
    elseif iname>2 && iname<6
        iInt = iname;
    else
        switch iname
            case 6
                iReg = 2;
            case 7
                isnr = 0.5;
            case 8
                isnr = 0.9;
            case 9
                iss=0;
            case 10
                iss = 0.25;
            case 11
                iss = 0.75;
            case 12
                iss = 1;
            case 13
                ilag = 1;
            case 14
                ifilt = 'e';
            case 15
                ifilt = 'c';
            case 16
                ifilt = 'cr';
        end
    end
    
    %%
    for iit= 1:nit
        tic
        iit
        clearvars -except iInt iReg isnr iss ilag ifilt iit nit iname DIRIN
        
        inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
            ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
        
        load([DIRIN inname '.mat'])
        
        %% define which pipelines run with this configuration
        if params.ip == 1
            params.pips = 1:20;
        elseif strcmp(params.ifilt,'c')|| strcmp(params.ifilt,'cr') %with champaign
            params.pips = [1:3 8];
        else % run only fixed, 90% and 99% pipelines
            params.pips = 1:9;
        end
        
        %%
        clear DIFFGC
        
        for ipip = params.pips
            
            %1 to 6: fixed
            %7: 90%
            %8: 99%
            %9: baseline
            %10: sumVox
            %11: 90% corrected
            %12: 99% corrected
            %13 to 20: equal to 1 to 8 but with zscoring
            %21: sum, then MIM
            %22: central voxel
            
            fprintf(['Testing pipeline ' num2str(ipip) '\n'])
            
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
                
                clear s_
                if strcmp(params.ifilt,'c') || strcmp(params.ifilt,'cr') %in case of the champaign filter
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
                        try
                            npcs(aroi) = min(find(varex> 0.99));
                            var_explained=0.99;
                        end
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
                    end
                    
                    if strcmp(params.ifilt,'c') || strcmp(params.ifilt,'cr') %champaign filter
                        %sort out rois with zero activity
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
                    %baseline and central voxel pipeline
                    signal_roi = cat(1,signal_roi,reshape(signal_source,[],l_epoch,n_trials));
                    npcs(aroi) = 3;
                    
                elseif ipip == 10
                    %sumVox
                    signal_roi{aroi} = reshape(signal_source,[],l_epoch,n_trials);
                    
                elseif ipip == 21
                    %sum of activity, then MIM
                    signal_roi = cat(1,signal_roi, squeeze(mean(reshape(signal_source,ndim, nvoxroi(aroi),l_epoch,n_trials),2)));
                    npcs(aroi) = 3;
                end
            end
            
            
            
            %% calculate connectivity scores
            
            if ipip ~= 10 && ipip ~= 11 && ipip ~= 12
                
                if strcmp(params.ifilt,'c') || strcmp(params.ifilt,'cr')
                    %no calculations for empty rois
                    npcs(empty_rois) = [];
                end
                
                % calculate indices
                clear inds PCA_inds
                [inds, PCA_inds] = fp_npcs2inds(npcs);
                
                % calculate MIM, MIC, TRGC and COH
                output = {'TRGC'};
                conn = data2sctrgcmim(signal_roi, fres, 50, 0,0, [], inds, output,0);
                
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
            
            %% save stuff
            
            MIM{ipip} = MIM_;
            MIC{ipip} = MIC_;
            
            if ipip ~= 10
                
                aCOH{ipip} = aCOH_;
                iCOH{ipip} = iCOH_;
                
                if ipip ~= 11 && ipip ~= 12 && ipip < 21
                    DIFFGC{ipip} = DIFFGC_;
                end
                
                if ipip ~= 9 && ipip < 21
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
                
            end
            
            
            
            %% Evaluate
            
            [mrr_mic(ipip), pr_mic(ipip),hk_mic(ipip),em1_mic(ipip),em2_mic(ipip),em3_mic(ipip)] ...
                = fp_mrr_hk(MIC_,iroi_seed,iroi_tar,1);
            [mrr_mim(ipip), pr_mim(ipip),hk_mim(ipip),em1_mim(ipip),em2_mim(ipip),em3_mim(ipip)] ...
                = fp_mrr_hk(MIM_,iroi_seed,iroi_tar,1);
            
            if ipip ~= 10
                [mrr_aCoh(ipip), pr_aCoh(ipip),hk_aCoh(ipip),em1_aCoh(ipip),em2_aCoh(ipip),em3_aCoh(ipip)] ...
                    = fp_mrr_hk(aCOH_,iroi_seed,iroi_tar,1);
                [mrr_iCoh(ipip), pr_iCoh(ipip),hk_iCoh(ipip),em1_iCoh(ipip),em2_iCoh(ipip),em3_iCoh(ipip)] ...
                    = fp_mrr_hk(iCOH_,iroi_seed,iroi_tar,1);
                
                if ipip ~= 11 && ipip ~= 12  && ipip < 21
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
                end
                
            end
            
            
            
            clear MIM_ MIC_ DIFFGC_ aCOH_ iCOH_
            
            
            
        end %pips
        
        %% Saving
        fprintf('Saving... \n')
        %save all
        outname = sprintf('%smim_%s.mat',DIROUT,params.logname);
        save(outname,'-v7.3')
        
        %save only evaluation parameters
        outname1 = sprintf('%smrr_%s.mat',DIROUT,params.logname);
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
end
        
        
        
        
        
        
        