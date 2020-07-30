function fp_mim_struct_sim(params,logname)

DIROUT = '/home/bbci/data/haufe/Franziska/data/';

if params.ip==7 && ~strcmp(params.ifilt,'l')
    load(sprintf('%s/mim_CS/filt/%d.mat',DIROUT,params.iit));
elseif params.ip==8 && ~params.ihemi==0
    load(sprintf('%s/mim_CS/hemi/%d.mat',DIROUT,params.iit));
else
    
    if params.ip==5 && ~params.iss==0
        load(sprintf('%s/mim_CS/noise_mix/%d.mat',DIROUT,params.iit));
    elseif params.ip==4 && ~params.isnr==0.1
        load(sprintf('%s/mim_CS/snr/%d.mat',DIROUT,params.iit));
    else
        
        fres = 40;
        n_trials = 200;
        
        %% signal generation
        
        % ROI labels
        % In D.sub_ind_roi, there are the randomly
        % selected voxels of each region
        fprintf('Getting atlas positions... \n')
        tic
        D = fp_get_Desikan(params.iReg);
        toc
        
        %signal generation
        fprintf('Signal generation... \n')
        tic
        [sig,brain_noise,sensor_noise,gt,L,iroi_seed, iroi_tar,D] = fp_generate_mim_signal(params, ...
            fres,n_trials, D,DIROUT);
        toc
        
        if params.ip==5 && params.iss==0
            outname = sprintf('%s/mim_CS/noise_mix/%d.mat',DIROUT,params.iit);
            save(outname,'-v7.3')
        elseif params.ip==4 && params.isnr==0.1
            outname = sprintf('%s/mim_CS/snr/%d.mat',DIROUT,params.iit);
            save(outname,'-v7.3')
        end
    end
    
    
    %combine noise sources
    for itrial = 1:n_trials
        noise = params.iss*brain_noise{itrial} + (1-params.iss)*sensor_noise{itrial};
        noise = noise ./ norm(noise, 'fro');
        %combine signal and noise
        signal_sensor1 = params.isnr*sig{itrial} + (1-params.isnr)*noise;
        signal_sensor(:,:,itrial) = signal_sensor1 ./ norm(signal_sensor1, 'fro');
    end
    
    
    %% get CS and filter A
    %parameters
    id_trials_1 = 1:n_trials;
    id_trials_2 = 1:n_trials;
    id_meg_chan = 1:size(signal_sensor,1);
    nmeg = numel(id_meg_chan);
    
    %cross spectrum
    fprintf('Calculating cross spectrum... \n')
    tic
    CS = fp_tsdata_to_cpsd(signal_sensor,fres,'WELCH',...
        [id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
    CS(:,:,1)=[];
    nfreq = size(CS,3);
    toc
    if params.ip==7 && strcmp(params.ifilt,'l')
        outname = sprintf('%s/mim_CS/filt/%d.mat',DIROUT,params.iit);
        save(outname,'-v7.3')
    elseif params.ip==8 && params.ihemi==0
        outname = sprintf('%s/mim_CS/hemi/%d.mat',DIROUT,params.iit);
        save(outname,'-v7.3')
    end
end

%% leadfield
L3 = L(:, D.ind_cortex, :);
for is=1:D.nvox
    clear L2
    L2 = L3(:,is,:);
    
    %remove radial orientation
    clear u s
    [u, s, v] = svd(squeeze(L2));
    L_backward(:,is,:) = u(:,:)*s(:,1:2);
end
ni = size(L_backward,3);

%construct source filter
if strcmp(params.ifilt,'e')
    A = squeeze(mkfilt_eloreta_v2(L_backward));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(params.ifilt,'d')
    
    A=zeros(nmeg,ni,D.nvox,nfreq);
    
    for ifrq = 1:nfreq
        cCS = CS(:,:,ifrq);
        lambda = mean(diag(real(cCS)))/100;
        
        CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
        
        for is=1:D.nvox %iterate across nodes
            Lloc=squeeze(L_backward(:,is,:));
            A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
        end
    end
    fqA = 1:nfreq; %This filter is frequency specific.
    nfqA = nfreq;
    
    
elseif strcmp(params.ifilt,'l')
    cCS = sum(CS,3);
    reg = 0.05*trace(cCS)/length(cCS);
    Cr = cCS + reg*eye(size(cCS,1));
    
    [~, A] = lcmv_meg(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;   
end

%% calculate MIM
%pca pipeline ('all' 8 pipelines + baseline)
[mic, mim, to_save] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,'all');

%% performance measures
fprintf('Performance measures... \n')
tic
[PERFORMANCE, BASELINE] = fp_get_performance(gt, mic, mim, params);
toc

%% save performance and baseline
fprintf('Saving... \n')
tic
outname = sprintf('%smim_%s.mat',DIROUT,logname);
save(outname,'-v7.3')
toc

