function fp_mim_struct_sim(params)

DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim3/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/mim_save3/';
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
    n_sensors = size(signal_sensor,1);
    
    %cross spectrum
    fprintf('Calculating cross spectrum... \n')
    CS = tsdata_to_cpsd_fast(signal_sensor,fres,'WELCH');
    CS(:,:,1)=[];
    nfreq = size(CS,3);
    
    if params.ip==1
       dir1 =  sprintf('%s/mim_CS/',DIROUT1);
       if ~exist(dir1); mkdir(dir1); end
       outname = sprintf('%s/mim_CS/%d.mat',DIROUT1,params.iit);
       save(outname,'-v7.3')
    end
    t.cs =toc;
end

%% leadfield
tic
L3 = L(:, D.ind_cortex, :);
L_backward = L3; 
ni = size(L_backward,3);

%construct source filter
if strcmp(params.ifilt,'e')
    reg_param = fp_eloreta_crossval(signal_sensor,L_backward,5);
    A = squeeze(mkfilt_eloreta_v2(L_backward,reg_param));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(params.ifilt,'d')
    
    A=zeros(n_sensors,ni,D.nvox,nfreq);
    
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
    cCS = sum(real(CS),3);
%     cCS = cov(signal_sensor(:,:)');
    reg = 0.05*trace(cCS)/length(cCS);
    Cr = cCS + reg*eye(size(cCS,1));
    
    [~, A] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;   
    
elseif strcmp(params.ifilt,'c')
    
%     sigma2_total = mean(diag(cov(signal_sensor(:, :)')));    
%     regu = sigma2_total*0.2;
%     sigu = regu*eye(n_sensors);
    sigu = cov(sensor_noise');
    tic
    L_backward = permute(L_backward,[1 3 2]);
    [~,~,w] = awsm_champ(signal_sensor(:, :), L_backward(:, :) ,...
        sigu, 200, 3, 2, 0);
    toc
    
    A = real(reshape(w',size(L_backward)));
%     A = permute(A,[1, 3, 2]);
    
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1; 
end

t.filter = toc;

%% calculate MIM

if 0%params.ip ==1
    
    %pca pipeline ('all' 8 pipelines + baseline)
    zs=1;
    [mic, mim, to_save, mean_icoh, mean_acoh,t] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,'all',zs,t);
    fprintf('pipelines calculated')

    
    %% corrected mim/mic
    fprintf('Correct MIM...\n')
    tic
    [mic,mim, to_save, mean_icoh,mean_acoh] = fp_correct_mim(A,signal_sensor, fqA, nfqA, D, params.ihemi, mic, mim, mean_icoh, mean_acoh, to_save,fres); 
    t.correctmims = toc;
    
    %% without ZS standardisation
    fprintf('Calculating fix, max and percent without ZS standardisation... \n')
    zs=0;
    for ii = 1:6
        [mic_fixed_zs{ii}, mim_fixed_zs{ii}, to_save_fixed_zs{ii},...
            mean_icoh_fixed_zs{ii},mean_acoh_fixed_zs{ii},t1] = ...
            fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,ii,zs,[]);
        t.zs0fixed{ii} = t1.mim;
    end
    [mic_max_zs, mim_max_zs, to_save_max_zs, mean_icoh_max_zs, mean_acoh_max_zs,t1] = ...
        fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,'max',zs,[]);
    t.zs0max = t1.mim;
    [mic_percent_zs, mim_percent_zs, to_save_percent_zs, mean_icoh_percent_zs,mean_acoh_percent_zs,t1] = ...
        fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,'percent',zs,[]);
    t.zs0percent = t1.mim;

    mic.fixed_zs0 = mic_fixed_zs;
    mim.fixed_zs0 = mim_fixed_zs; 
    to_save.fixed_zs0 = to_save_fixed_zs;
    mean_icoh.fixed_zs0 = mean_icoh_fixed_zs; 
    mean_acoh.fixed_zs0 = mean_acoh_fixed_zs; 

    mic.max_zs0 = mic_max_zs; 
    mim.max_zs0 = mim_max_zs; 
    to_save.max_zs0 = to_save_max_zs;
    mean_icoh.max_zs0 = mean_icoh_max_zs; 
    mean_acoh.max_zs0 = mean_acoh_max_zs; 

    mic.percent_zs0 = mic_percent_zs; 
    mim.percent_zs0 = mim_percent_zs; 
    to_save.percent_zs0 = to_save_percent_zs;
    mean_icoh.percent_zs0 = mean_icoh_percent_zs;
    mean_acoh.percent_zs0 = mean_acoh_percent_zs;

    clear mic_max_zs mim_max_zs to_save_max_zs mean_icoh_max_zs mean_acoh_max_zs mic_percent_zs ...
        mim_percent_zs to_save_percent_zs mean_icoh_percent_zs mean_acoh_percent_zs mic_fixed_zs ...
        mim_fixed_zs to_save_fixed_zs mean_icoh_fixed_zs mean_acoh_fixed_zs
else
    %reduced pca pipelines: zs0 of all fixed pips and of 90 and 99 %,
    %baseline 
    zs=0;
    [mic, mim, to_save, mean_icoh, mean_acoh,t] = fp_get_mim_reduced(A,CS,fqA,nfqA, D,params.ihemi,'all',zs,t);
    fprintf('pipelines calculated')
end

to_save.t = t;

%% save performance and baseline
fprintf('Saving... \n')
tic
outname = sprintf('%smim_%s.mat',DIROUT,params.logname);
save(outname,'-v7.3')
t.saving = toc;

