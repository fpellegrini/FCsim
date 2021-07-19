function fp_add_corrected_pips

fp_addpath

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim4/';
DIRLOG ='/home/bbci/data/haufe/Franziska/log/mim_sim4/corrected/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

iit = str2num(getenv('SGE_TASK_ID'));
%%

%default paramenters
iInt = 2;
iReg=1;
isnr=0.7;
iss = 0.5;
ilag=2;
ifilt='l';
logname = sprintf('iter%d',iit);


%%

if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
    eval(sprintf('!touch %s%s_work',DIRLOG,logname))
    fprintf('Working on %s. \n',logname)
    
    tic
    
    iit
    
    inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
        ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
    
    load([DIRIN inname '.mat'])
    
    for ipip = 11:12
        
        fprintf(['Testing pipeline ' num2str(ipip) '\n'])
        
        tic
        %% PCA
        
        clear npcs
        signal_roi = []
        
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
            s_ = logical(ones(size(signal_source,1),1));
            
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
                try %might be empty in case of the cr filter
                    npcs(aroi) = min(find(varex> 0.99));
                catch
                    npcs(aroi) = 0;
                end
                var_explained=0.99;
            end
            
            %bring signal_roi to the shape of npcs x l_epoch x n_trials
            signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
        end
        
        
        %% calculate connectivity scores
        
        % calculate indices
        clear inds PCA_inds
        [inds, PCA_inds] = fp_npcs2inds(npcs);
        
        %calculate normalized values of MIM and MIC
        [MIM_, MIC_, iCOH_, aCOH_] = fp_correct_mim(signal_roi,inds,fres,D,filt,PCA_inds);
        
        
        t.pips(ipip) = toc;
        
        %% save stuff
        
        MIM{ipip} = MIM_;
        MIC{ipip} = MIC_;
        
        aCOH{ipip} = aCOH_;
        iCOH{ipip} = iCOH_;
        
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
    
    eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
    
end

