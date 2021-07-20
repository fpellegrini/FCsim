function fp_data2mim_sim_dics(params)

% define folders for saving results
DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim4/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/mim_save4/';
if ~exist(DIROUT1);mkdir(DIROUT1); end

%define which pipelines run with this configuration
params.pips = 1:9;

%% load data

%reload data from ip1 to keep them constant and only vary the source
%localization
params_save = params;
load(sprintf('%s/mim_sensorsig/%d.mat',DIROUT1,params.iit));
params = params_save;
clear params_save

%% leadfield and inverse solution

tic

%select only voxels that belong to any roi
L_backward = L(:, D.ind_cortex, :);
ndim = size(L_backward,3);

%CS
CSpara = [];
CSpara.subave = 0;
CSpara.mywindow = hanning(l_epoch)./sqrt(hanning(l_epoch)'*hanning(l_epoch));
CS_sensor = data2cs_event(signal_sensor(:, :)', l_epoch, l_epoch, l_epoch, fres+1, CSpara);

%construct source filter
A1 = [];
for ifreq = 1:fres+1
    cCS = real(CS_sensor(:, :, ifreq));
    reg = 0.05*trace(cCS)/length(cCS);
    Cr = cCS + reg*eye(n_sensors);
    
    [~, A1(:, :, :, ifreq)] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));    
end
A1 = permute(A1,[1, 3, 2, 4]);

t.filter = toc;

%% Loop over different pipelines

errorpipeline = [];
for ipip = params.pips %most successful: ipip 1 to 3
    
    %1 to 6: fixed
    %7: 90%
    %8: 99%
    %9: baseline
    
    fprintf(['Testing pipeline ' num2str(ipip) '\n'])
    try
        tic
        %% PCA
        
        clear npcs P A2 V CS_roi
        
        %loop over regions
        for aroi = 1:D.nroi
            
            clear A_ CS_source
            
            %A_ is the lcmv filter at aroi
            A_ = A1(:, :,D.ind_roi_cortex{aroi},:);
            
            if ipip == 9 %baseline: pre-select voxels of true activity
                A_ = A_(:,:,D.sub_ind_roi_region{aroi},:);
            end
            
            %number of voxels at the current roi
            nvoxroi(aroi) = size(A_,3);
            A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi),fres+1]);
            
            %project sensor signal to voxels at the current roi (aroi)
            for ifq = 1:fres+1
                CS_source(:,:,ifq) = squeeze(A2{aroi}(:,:,ifq))' * squeeze(CS_sensor(:,:,ifq))...
                    * squeeze(A2{aroi}(:,:,ifq));
            end
            
            
            if ~(ipip == 9) && ~(ipip == 10) && ipip < 21
                
                %PCA
                clear CSs v v5 in V_ D_
                CSs = squeeze(sum(real(CS_source),3)); %covariance
                [V_, D_] = eig(CSs);
                [D_, in] = sort(real(diag(D_)), 'descend');
                V{aroi} = V_(:,in);
                vx_ = cumsum(D_)./sum(D_);
                
                % variance explained
                invx = 1:min(length(vx_), n_sensors);
                varex = vx_(invx);
                
                %save npcs and variance explained
                if ismember(ipip,[7])
                    %npcs are selected in a way that 90% of the variance is preserved
                    npcs(aroi) = min(find(varex> 0.9));
                    var_explained=0.9;
                elseif ismember(ipip,[8])
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
                end
                
                V{aroi} = real(V{aroi}(:, 1:npcs(aroi)));
                
                for ifq = 1:fres+1
                    P(:,:,aroi,ifq) = squeeze(A2{aroi}(:,:,ifq)) * V{aroi}; 
                end  
                
            elseif ipip == 9
                %baseline
                P(:,:,aroi,:) = A2{aroi};
                npcs(aroi) = nvoxroi(aroi)*ndim;
            end
           
        end
        
        CS_roi = [];
        for ifreq = 1:fres+1
            CS_roi(:, :, ifreq) = reshape(P(:, :, :, ifreq), n_sensors, [])' ...
                * CS_sensor(:, :, ifreq) * reshape(P(:, :, :, ifreq), n_sensors, []);
        end
        
        
        
        %% calculate connectivity scores        
        
        % calculate indices
        clear inds PCA_inds
        [inds, PCA_inds] = fp_npcs2inds(npcs);
        
        output = {'MIM','MIC','TRGC','COH'};
        
        %calculation of output of data2sctrgcmim based on CS
        conn = fp_cs2sctrgcmim(CS_roi, fres, 30, inds, output);
        
        % extract measures out of the conn struct
        [MIM_, MIC_, DIFFGC_, iCOH_, aCOH_] = fp_unwrap_conn(conn,numel(npcs),filt,PCA_inds);        
   
        t.pips(ipip) = toc;
        
        %% save stuff
        
        MIM{ipip} = MIM_;
        MIC{ipip} = MIC_;
        
        if ipip ~= 10
            
            aCOH{ipip} = aCOH_;
            iCOH{ipip} = iCOH_;
            
            if ipip ~= 11 && ipip ~= 12
                DIFFGC{ipip} = DIFFGC_;
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
                    if ~strcmp(params.ifilt,'c') || ~strcmp(params.ifilt,'cr')
                        corr_voxnpcs(ipip) = corr(nvoxroi', npcs');
                    end
                end
            end
            
        end
        
        
        
        %% Evaluate
        
        [mrr_mic(ipip), pr_mic(ipip),hk_mic(ipip),em3_mic(ipip)] ...
            = fp_mrr_hk(MIC_,iroi_seed,iroi_tar,1);
        [mrr_mim(ipip), pr_mim(ipip),hk_mim(ipip),em3_mim(ipip)] ...
            = fp_mrr_hk(MIM_,iroi_seed,iroi_tar,1);
        
        if ipip ~= 10
            [mrr_aCoh(ipip), pr_aCoh(ipip),hk_aCoh(ipip),em3_aCoh(ipip)] ...
                = fp_mrr_hk(aCOH_,iroi_seed,iroi_tar,1);
            [mrr_iCoh(ipip), pr_iCoh(ipip),hk_iCoh(ipip),em3_iCoh(ipip)] ...
                = fp_mrr_hk(iCOH_,iroi_seed,iroi_tar,1);
            
            if ipip ~= 11 && ipip ~= 12
                %absolute value of gc and only triu is considered. Metric neglects
                %the direction of the interaction
                [mrr_absgc(ipip), pr_absgc(ipip),hk_absgc(ipip),em3_absgc(ipip)] ...
                    = fp_mrr_hk(abs(DIFFGC_),iroi_seed,iroi_tar,1);
                
                %only positive part of gc is submitted and the whole matrix is
                %considered. Metric that is strongly influenced by the direction of
                %the effect
                clear pos_diffgc
                pos_diffgc = DIFFGC_;
                pos_diffgc(pos_diffgc< 0) = 0;
                [mrr_posgc(ipip), pr_posgc(ipip),hk_posgc(ipip), em3_posgc(ipip)] ...
                    = fp_mrr_hk(pos_diffgc,iroi_seed,iroi_tar,0);
            end
            
        end
        
        
        
        clear MIM_ MIC_ DIFFGC_ aCOH_ iCOH_
        
    catch
        errorpipeline = [errorpipeline ipip];       
    end
    
end %pips

%% Saving
fprintf('Saving... \n')
%save all
outname = sprintf('%smim_%s.mat',DIROUT,params.logname);
save(outname,'-v7.3')

%save only evaluation parameters
outname1 = sprintf('%smrr_%s.mat',DIROUT,params.logname);
save(outname1,...
    'mrr_mic','pr_mic','hk_mic','em3_mic',...
    'mrr_mim','pr_mim','hk_mim','em3_mim',...
    'mrr_aCoh','pr_aCoh','hk_aCoh','em3_aCoh',...
    'mrr_iCoh','pr_iCoh','hk_iCoh','em3_iCoh',...
    'mrr_absgc','pr_absgc','hk_absgc','em3_absgc',...
    'mrr_posgc','pr_posgc','hk_posgc','em3_posgc',...
    '-v7.3')

