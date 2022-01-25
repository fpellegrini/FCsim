function fp_addpip

fp_addpath

DIRLOG ='/home/bbci/data/haufe/Franziska/log/mim_sim5/mean+FC/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim5/';
DIROUT_new = '/home/bbci/data/haufe/Franziska/data/mim_sim5/addmean+FC/';
if ~exist(DIROUT_new); mkdir(DIROUT_new); end


%%

for iit = 1:100
    
    clearvars -except iit DIRIN DIRLOG
    
    %default paramenters
    iInt = 2;
    iReg=1;
    isnr=0.3;
    iss = 0.5;
    ilag=2;
    ifilt='l';
    dimred = 'p'; %pca
    
    inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_iter%d'...
        ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,dimred,iit);
    
    load([DIRIN inname '.mat'])
    
    if ~exist(sprintf('%s%s_work',DIRLOG,params.logname)) & ~exist(sprintf('%s%s_done',DIRLOG,params.logname))
        eval(sprintf('!touch %s%s_work',DIRLOG,params.logname))
        fprintf('Working on %s. \n',params.logname)
        
        for ipip = 21
            
            %1 to 6: fixed
            %7: 90%
            %8: 99%
            %9: baseline
            %10: sumVox
            %11: 90% corrected
            %12: 99% corrected
            %13 to 20: equal to 1 to 8 but WITHOUT zscoring
            %21: mean, then MIM
            %22: central voxel
            
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
                    
                    if ipip == 9 %baseline: pre-select voxels of true activity
                        A_ = A(:, :,D.ind_roi_cortex{aroi},:);
                        A_ = A_(:,:,D.sub_ind_roi_region{aroi},:);
                        
                    elseif ipip == 22 %pick central voxel of each region
                        A_ = A(:, :,D.ind_roi_cortex{aroi},:);
                        A_ = A_(:,:,D.ctr_ind_roi_region{aroi},:);
                        
                    else
                        A_ = A(:, :,D.ind_roi_cortex{aroi},:);
                    end
                    
                    %number of voxels at the current roi
                    nvoxroi(aroi) = size(A_,3);
                    A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
                    
                    %project sensor signal to voxels at the current roi (aroi)
                    signal_source = A2{aroi}' * signal_sensor(:,:);
                    
                    if ipip < 13 && ipip > 20 %zscoring all pipelines but 13 to 20
                        signal_source = zscore(signal_source);
                    end
                    
                    clear s_
                    if strcmp(params.ifilt(1),'c') %in case of the champaign filter
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
                            %do SSD
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
                        %sumVox
                        signal_roi{aroi} = reshape(signal_source,[],l_epoch,n_trials);
                        
                    elseif ipip == 21
                        %mean of activity, then MIM
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
                    
                    % calculate MIM, MIC, TRGC and COH
                    %             if ipip > 20
                    %                 output = {'MIM','MIC','COH'};
                    %             else
                    output = {'MIM','MIC','GC','TRGC','COH'};
                    %             end
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
                    
                    if ipip ~= 11 && ipip ~= 12
                        DIFFGC{ipip} = DIFFGC_;
                        GC{ipip} = GC_;
                    end
                    
                    if ipip ~= 9
                        try
                            to_save{ipip}.npcs = npcs;
                            to_save{ipip}.varex = var_explained;
                            
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
                
            catch
                errorpipeline = [errorpipeline ipip];
            end
            
        end %pips
        
        %% Saving
        fprintf('Saving... \n')
        %save all
        outname = sprintf('%smim_%s.mat',DIROUT_new,params.logname);
        save(outname,'-v7.3')
        
        %save only evaluation parameters
        outname1 = sprintf('%spr_%s.mat',DIROUT_new,params.logname);
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
        
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,params.logname,DIRLOG,params.logname))
    end
end