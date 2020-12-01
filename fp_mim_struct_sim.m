function fp_mim_struct_sim(params)

DIROUT = '/home/bbci/data/haufe/Franziska/data/mim_sim1/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/mim_save1/';
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
        
        % ROI labels
        % In D.sub_ind_roi, there are the randomly
        % selected voxels of each region
        fprintf('Getting atlas positions... \n')
        tic
        D = fp_get_Desikan(params.iReg);
        toc
        
        %signal generation
        fprintf('Signal generation... \n')
        [sig,brain_noise,sensor_noise,gt,L,iroi_seed, iroi_tar,D, fres, n_trials] = fp_generate_mim_signal(params, ...
            D,DIROUT1);
     
        if params.ip==1
            dir1 =  sprintf('%s/mim_sig/',DIROUT1);
            if ~exist(dir1); mkdir(dir1); end
            outname = sprintf('%s/mim_sig/%d.mat',DIROUT1,params.iit);
            save(outname,'-v7.3')
        end
    end
       
    %combine noise sources
    noise = params.iss*brain_noise + (1-params.iss)*sensor_noise;
    noise = noise ./ norm(noise, 'fro');
    %combine signal and noise
    signal_sensor1 = params.isnr*sig + (1-params.isnr)*noise;
    signal_sensor = signal_sensor1 ./ norm(signal_sensor1, 'fro');
    
    %reshape
    signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
    
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
    
    if params.ip==1
       dir1 =  sprintf('%s/mim_CS/',DIROUT1);
       if ~exist(dir1); mkdir(dir1); end
       outname = sprintf('%s/mim_CS/%d.mat',DIROUT1,params.iit);
       save(outname,'-v7.3')
    end
end

%% leadfield
L3 = L(:, D.ind_cortex, :);
L_backward = L3; 
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
    
    [~, A] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;   
end

%% calculate MIM
%pca pipeline ('all' 8 pipelines + baseline)
zs=1;
[mic, mim, to_save, mean_coh] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,'all',zs);
fprintf('pipelines calculated')
d=whos; sum([d.bytes])/1000^3

%% corrected mim/mic

[mim,mic,mean_coh,to_save] = fp_correct_mim(A,signal_sensor, fqA, nfqA, D, params.ihemi, mic, mim, mean_coh, to_save); 

%% without ZS standardisation
zs=0;
for ii = 1:5
    [mic_fixed_zs{ii}, mim_fixed_zs{ii}, to_save_fixed_zs{ii}, mean_coh_fixed_zs{ii}] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,ii,zs);
end
[mic_max_zs, mim_max_zs, to_save_max_zs, mean_coh_max_zs] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,'max',zs);
[mic_percent_zs, mim_percent_zs, to_save_percent_zs, mean_coh_percent_zs] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,'percent',zs);

mic.fixed_zs0 = mic_fixed_zs;
mim.fixed_zs0 = mim_fixed_zs; 
to_save.fixed_zs0 = to_save_fixed_zs;
mean_coh.fixed_zs0 = mean_coh_fixed_zs; 

mic.max_zs0 = mic_max_zs; 
mim.max_zs0 = mim_max_zs; 
to_save.max_zs0 = to_save_max_zs;
mean_coh.max_zs0 = mean_coh_max_zs; 

mic.percent_zs0 = mic_percent_zs; 
mim.percent_zs0 = mim_percent_zs; 
to_save.percent_zs0 = to_save_percent_zs;
mean_coh.percent_zs0 = mean_coh_percent_zs;

clear mic_max_zs mim_max_zs to_save_max_zs mean_coh_max_zs mic_percent_zs ...
    mim_percent_zs to_save_percent_zs mean_coh_percent_zs mic_fixed_zs ...
    mim_fixed_zs to_save_fixed_zs mean_coh_fixed_zs

%% performance measures
fprintf('Performance measures... \n')
tic
[PERFORMANCE, BASELINE] = fp_get_performance(gt, mic, mim, mean_coh);
toc
d=whos; sum([d.bytes])/1000^3

%% save performance and baseline
fprintf('Saving... \n')
tic
outname = sprintf('%smim_%s.mat',DIROUT,params.logname);
save(outname,'-v7.3')
toc

