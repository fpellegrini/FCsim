function fp_data2mim_sim(params)

% define folders for saving results 
DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim3/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/mim_save3/';
if ~exist(DIROUT1);mkdir(DIROUT1); end


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
        
        if params.ip==1
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
    
    signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
    signal_sensor = signal_sensor / norm(signal_sensor, 'fro');
    
    %reshape
    signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
    [n_sensors, l_epoch, n_trials] = size(signal_sensor);
    
    t.snr = toc;
    
    if params.ip==1
        dir1 =  sprintf('%s/mim_sensorsig/',DIROUT1);
        if ~exist(dir1); mkdir(dir1); end
        outname = sprintf('%s/mim_sensorsig/%d.mat',DIROUT1,params.iit);
        save(outname,'-v7.3')
    end
    
end

%% leadfield and inverse solution

tic
L_backward = L(:, D.ind_cortex, :);
ni = size(L_backward,3);

%construct source filter
if strcmp(params.ifilt,'e')
    reg_param = fp_eloreta_crossval(signal_sensor,L_backward,5);
    A = squeeze(mkfilt_eloreta_v2(L_backward,reg_param));
    A = permute(A,[1, 3, 2]);
    
elseif strcmp(params.ifilt,'l')
    cCS = cov(signal_sensor(:,:)');
    reg = 0.05*trace(cCS)/length(cCS);
    Cr = cCS + reg*eye(size(cCS,1));
    
    [~, A] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
    A = permute(A,[1, 3, 2]);
    
elseif strcmp(params.ifilt,'c')
    
    %     sigma2_total = mean(diag(cov(signal_sensor(:, :)')));
    %     regu = sigma2_total*0.2;
    %     sigu = regu*eye(n_sensors);
    sigu = cov(sensor_noise');
    L_backward = permute(L_backward,[1 3 2]);
    [~,~,w] = awsm_champ(signal_sensor(:, :), L_backward(:, :) ,...
        sigu, 200, 3, 2, 0);
    
    A = real(reshape(w',size(L_backward)));
end

t.filter = toc;


%% Loop over different pipelines

for ipip = 1:20
    
    %1 to 6: fixed
    %7: 90%
    %8: 99%
    %9: baseline
    %10: sumVox
    %11: 90% corrected
    %12: 99% corrected
    %13 to 20: equal to 1 to 8 but with zscoring
    
    clear MIM_ MIC_ iCOH_ aCOH_ 
    tic
    %% PCA
    
    clear npcs
    signal_roi = [];
    
    %loop over regions
    for aroi = 1:D.nroi
        
        clear A_ signal_source
        
        %A_ is the lcmv filter at aroi
        A_ = A(:, :,D.ind_roi_cortex{aroi},:);
        
        if ipip == 9
            %baseline: pre-select voxels of true activity
            A_ = A_(:,:,D.sub_ind_roi_region{aroi});
        end
        
        %number of voxels at the current roi
        nvoxroi(aroi) = size(A_,3);
        
        A2{aroi} = reshape(A_, [n_sensors, ni*nvoxroi(aroi)]);
        
        %project sensor signal to voxels at the current roi (aroi)
        signal_source = A2{aroi}' * signal_sensor(:,:);
        
        %zscoring
        if ipip > 12
            signal_source = zscore(signal_source)*sqrt(mean(var(signal_source)));
        end
        
        if ~(ipip == 9) && ~(ipip == 10)
            %do PCA
            clear signal_roi_ S
            [signal_roi_,S,~] = svd(double(signal_source)','econ');
            
            % variance explained
            vx_ = cumsum(diag(S).^2)./sum(diag(S).^2);
            invx = 1:min(length(vx_), n_sensors);
            varex = vx_(invx);
            
            if ismember(ipip,[7 11])
                npcs(aroi) = min(find(varex> 0.9));
                var_explained=0.9;
            elseif ismember(ipip,[8 12])
                npcs(aroi) = min(find(varex> 0.99));
                var_explained=0.90;
            elseif ipip <= 6
                npcs(aroi) = ipip;
                var_explained(aroi) = varex(npcs(aroi));
            elseif ipip >= 13
                npcs(aroi) = ipip-12;
                var_explained(aroi) = varex(npcs(aroi));
            end
            
            %bring signal_roi to the shape of npcs x l_epoch x n_trials
            signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi))* S(1:npcs, 1:npcs))',[],l_epoch,n_trials));
            
        elseif ipip == 9
            %baseline
            signal_roi = cat(1,signal_roi,reshape(signal_source,[],l_epoch,n_trials));
            npcs(aroi) = 3; %%%%%%%%%%%%%
            var_explained = [];
            
        elseif ipip == 10
            %sumVox
            signal_roi{aroi} = reshape(signal_source,[],l_epoch,n_trials);
        end
    end
    
    
    
    %% calculate MIM/MIC
   
    
    % Calculate mim, mic, iCoh, aCoh and trgc for all pipelines besides sumVox.
    if ipip ~= 10
        
        %% calculate indices
        clear beg_inds end_inds PCA_inds conn output
        
        beg_inds = cumsum([1 npcs(1:end-1)]);
        end_inds = cumsum([npcs]);
        
        for iroi = 1:numel(npcs)
            PCA_inds{iroi} = beg_inds(iroi):end_inds(iroi);
        end
        
        inds = {}; ninds = 0;
        for iroi = 1:numel(npcs)
            for jroi = (iroi+1):numel(npcs)
                inds{ninds+1} = {PCA_inds{iroi}, PCA_inds{jroi}};
                ninds = ninds + 1;
            end
        end
        
        
        output = {'MIM','MIC','TRGC','COH'};
        conn = data2sctrgcmim(signal_roi, fres, 20, 0,0, [], inds, output);
        
        % extract measures out of the conn struct
        iinds = 0;
        for iroi = 1:D.nroi
            for jroi = (iroi+1):D.nroi
                iinds = iinds + 1;
                MIM_(iroi, jroi) = mean(conn.MIM(filt.band_inds, iinds));
                MIC_(iroi, jroi) = mean(conn.MIC(filt.band_inds, iinds));
                DIFFGC_(iroi,jroi) = mean(squeeze(conn.TRGC(filt.band_inds,iinds,1) - conn.TRGC(filt.band_inds,iinds,2)));
            end
        end
        
        clear iCOH aCOH
        for iroi = 1:D.nroi
            for jroi = 1:D.nroi
                iCOH_(iroi, jroi) = mean(mean(mean(abs(imag(conn.COH(filt.band_inds, PCA_inds{iroi}, PCA_inds{jroi}))),1), 2), 3);
                aCOH_(iroi, jroi) = mean(mean(mean(abs(conn.COH(filt.band_inds, PCA_inds{iroi}, PCA_inds{jroi})),1), 2), 3);
            end
        end
        
    elseif ipip == 10
        % In sumVox,trgc would take too long
        output = {'MIM','MIC'};
        
        % For memory purposes, we need to calculate this region by region
        for oroi = 1:D.nroi-1
            for uroi = oroi+1:D.nroi
                clear data npcs beg_inds end_inds PCA_inds mic mim conn 
                
                npcs = repmat(ni,1,nvoxroi(oroi)+nvoxroi(uroi));
                data = cat(1,signal_roi{oroi}, signal_roi{uroi});
                
                beg_inds = cumsum([1 npcs(1:end-1)]);
                end_inds = cumsum([npcs]);
                
                for iroi = 1:numel(npcs)
                    PCA_inds{iroi} = beg_inds(iroi):end_inds(iroi);
                end
                
                inds = {}; ninds = 0;
                for iroi = 1:numel(npcs)
                    for jroi = (iroi+1):numel(npcs)
                        inds{ninds+1} = {PCA_inds{iroi}, PCA_inds{jroi}};
                        ninds = ninds + 1;
                    end
                end
                
                conn = data2sctrgcmim(data, fres, 20, 0,0, [], inds, output);
                
                % extract measures out of the conn struct
                iinds = 0;
                for ivox = 1:nvoxroi(oroi)+nvoxroi(uroi)-1
                    for jvox = (ivox+1):nvoxroi(oroi)+nvoxroi(uroi)
                        iinds = iinds + 1;
                        mim(ivox, jvox,:) = conn.MIM(:, iinds)';
                        mic(ivox, jvox,:) = conn.MIC(:, iinds)';
                    end
                end
                
                MIM_(oroi,uroi) = squeeze(mean(mean(mean(mim(1:nvoxroi(oroi),nvoxroi(oroi):end,filt.band_inds),1),2),3));
                MIC_(oroi,uroi) = squeeze(mean(mean(mean(mic(1:nvoxroi(oroi),nvoxoig(oroi):end,filt.band_inds),1),2),3));
                aCOH_=[];
                iCOH_=[];
                
            end
        end
        npcs = [];
        
        
    end %ipip == 10
    
    t.pips(ipip) = toc;
    
    %% save a few things 
    
    MIM{ipip} = MIM_;
    MIC{ipip} = MIC_;
    
    if ipip ~= 10
        
        DIFFGC{ipip} = DIFFGC_;
        aCOH{ipip} = aCOH_;
        iCOH{ipip} = iCOH_;
        
        if ipip ~= 9
            to_save{ipip}.npcs = npcs;
            to_save{ipip}.varex = var_explained;

            nvoxroi_all = nvoxroi'*nvoxroi;
            nvoxroi_all = nvoxroi_all(:);
            to_save{ipip}.corr_voxmim = corr(nvoxroi_all,MIM_(:));
            to_save{ipip}.corr_voxmic = corr(nvoxroi_all ,MIC_(:));
            to_save{ipip}.corr_voxicoh = corr(nvoxroi_all,iCOH_(:));
            to_save{ipip}.corr_voxacoh = corr(nvoxroi_all,aCOH_(:));
            to_save{ipip}.corr_voxnpcs = corr(nvoxroi', npcs');
        end
        
    end
    
    clear MIM_ MIC_ DIFFGC_ aCOH_ iCOH_ 


    
    
    
end %pips

fprintf('Saving... \n')
outname = sprintf('%smim_%s.mat',DIROUT,params.logname);
save(outname,'-v7.3')







%%

% correct

