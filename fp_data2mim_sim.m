function fp_data2mim_sim(params)

DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim4/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/mim_save4/';
if ~exist(DIROUT1);mkdir(DIROUT1); end

if params.ip==7 || params.ip==8
    params_save = params;
    load(sprintf('%s/mim_CS/%d.mat',DIROUT1,params.iit));
    params = params_save;
    clear params_save
else
    
    if params.ip==5 || params.ip ==4
        params_save = params;
        load(sprintf('%s/mim_sig/%d.mat',DIROUT1,params.iit));
        params= params_save;
        clear params_save
    else
        
        %% signal generation
        tic
        % ROI labels
        % In D.sub_ind_roi, there are the randomly
        % selected voxels of each region
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
    
    signal_sensor = filtfilt(filt.bhigh, filt.ahigh, signal_sensor1);
    signal_sensor = signal_sensor / norm(signal_sensor, 'fro');
    
    %reshape
    signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
    
    t.snr = toc;
    
    %% get CS and filter A
    tic
    %parameters
    [n_sensors, l_epoch, n_trials] = size(signal_sensor);
    
    if params.ip==1
        dir1 =  sprintf('%s/mim_CS/',DIROUT1);
        if ~exist(dir1); mkdir(dir1); end
        outname = sprintf('%s/mim_CS/%d.mat',DIROUT1,params.iit);
        save(outname,'-v7.3')
    end
    t.cs =toc;
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
    tic
    [~,~,w] = awsm_champ(signal_sensor(:, :), L_backward(:, :) ,...
        sigu, 200, 3, 2, 0);
    toc
    
    A = reshape(w',size(L_backward));
    A = permute(A,[1, 3, 2]);
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
    
    %%
    if ~(ipip==10)
        
        %% PCA
        
        %loop over regions
        for aroi = 1:D.nroi
            
            clear A_ signal_source
            
            %A_ is the lcmv filter at aroi
            A_ = A(:, :,D.ind_roi_cortex{aroi},:);
            
            if ipip == 9 
                %pre-select voxels of true activity 
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
            
            if ~(ipip == 9)
                %do PCA
                clear signal_roi_ S_
                [signal_roi_,S,~] = svd(double(signal_source)','econ');

                % variance explained
                vx_ = cumsum(diag(S).^2)./sum(diag(S).^2);
                invx = 1:min(length(vx_), n_sensors);
                varex = vx_(invx);

                if ismember(ipip,[7 11])
                    npcs(aroi) = min(find(varex> 0.9));
                elseif ismember(ipip,[8 12])
                    npcs(aroi) = min(find(varex> 0.99));
                elseif ipip <= 6
                    npcs(aroi) = ipip;
                elseif ipip >= 13
                    npcs(aroi) = ipip-12;
                end

                %bring signal_roi to the shape of npcs x l_epoch x n_trials
                signal_roi{aroi} = reshape(signal_roi_(:,1:npcs(aroi))',[],l_epoch,n_trials);
            else
                signal_roi{aroi} = reshape(signal_source,[],l_epoch,n_trials);
            end
        end
        
        
        %% calculate MIM/MIC
        
        %loop over all roi combinations
        tic
        for oroi = 1:nroi-1
            for uroi = oroi+1:nroi
                
                fprintf(['Calculating regions ' num2str(oroi) ' and ' num2str(uroi) '.\n'])
                
                clear data
                data = cat(1,signal_roi{oroi},signal_roi{uroi});
                
                clear mim mic trgc
                inds = fp_npcs2inds([npcs(oroi) npcs(uroi)]);
                [trgc,~,mim,mic, ~, inds] = data2sctrgcmim(data, fres, [], 0, 0, [], inds);
                
                %so far, mim{1} is equal to mim{2}
                MIM{ipip}([oroi, uroi],[oroi, uroi],:)=mim{1};
                MIC{ipip}([oroi, uroi],[oroi, uroi],:)=mic{1};
                DIFFGC{ipip}(oroi,uroi,:) = squeeze(trgc(:,1) - trgc(:,2));
                
                %mean iCoh and mean aCoh
                clear CS COH
                CS = tsdata_to_cpsd_fast(data, fres, 'WELCH', l_epoch);
                
                for ifreq = 1: fres+1
                    clear pow
                    pow = real(diag(CS(:,:,ifreq)));
                    COH(:,:,ifreq) = CS(:,:,ifreq)./ sqrt(pow*pow');
                end
                
                MEAN_ACOH{ipip}([oroi, uroi],[oroi, uroi],:)= mean(mean(abs(COH(inds{1,1}{1,1},inds{1,1}{1,2},:)),1),2);
                MEAN_ICOH{ipip}([oroi, uroi],[oroi, uroi],:)= mean(mean(abs(imag(COH(inds{1,1}{1,1},inds{1,1}{1,2},:))),1),2);
                
            end
        end
        toc
        
    elseif ipip == 10
        
        for aroi = 1:D.nroi
            clear A_ signal_source
            
            %A_ is the lcmv filter at aroi
            A_ = A(:, :,D.ind_roi_cortex{aroi});

            %number of voxels at the current roi
            nvoxroi(aroi) = size(A_,3);
            
            A2{aroi} = reshape(A_, [n_sensors, ni*nvoxroi(aroi)]);
            
            %project sensor signal to voxels at the current roi (aroi)
            signal_source{aroi} = A2{aroi}' * signal_sensor(:,:);
            
            npcs(aroi) = ni;
        end
        
       
                    
                end
                
                if ipip == 10
                    
                    %MIC and MIM
                    npcs = repmat(ndim,size(Cohroi,1)/(ndim),1);
                    [mic_v,mim_v] =  fp_mim(Cohroi,npcs);
                    
                    %sum up mim/ mic within rois
                    mic2(oroi,uroi,:) = squeeze(mean(mean(mic_v(1:nvoxreg1,nvoxreg1+1:end,:),1),2));
                    mic2(uroi,oroi,:) = mic2(oroi,uroi,:);
                    mic2(oroi,oroi,:) = squeeze(mean(mean(mic_v(1:nvoxreg1,1:nvoxreg1,:),1),2));
                    mim2(oroi,uroi,:) = squeeze(mean(mean(mim_v(1:nvoxreg1,nvoxreg1+1:end,:),1),2));
                    mim2(uroi,oroi,:) = mim2(oroi,uroi,:);
                    mim2(oroi,oroi,:) = squeeze(mean(mean(mim_v(1:nvoxreg1,1:nvoxreg1,:),1),2));
                    mean_icoh =[];
                    
                end
            end
        end
        
    end
end












%%

% sumVox,baseline, correct
% t, to_save mit corrs etc

