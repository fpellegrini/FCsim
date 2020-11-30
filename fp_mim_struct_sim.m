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
    
    %% reshape here! 
    signal_sensor2 = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
    
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
    
%     if params.ip==1
%        dir1 =  sprintf('%s/mim_CS/',DIROUT1);
%        if ~exist(dir1); mkdir(dir1); end
%        outname = sprintf('%s/mim_CS/%d.mat',DIROUT1,params.iit);
%        save(outname,'-v7.3')
%     end
end

%% leadfield
L3 = L(:, D.ind_cortex, :);
L_backward = L3; 
% for is=1:D.nvox
%     clear L2
%     L2 = L3(:,is,:);
%     
%     %remove radial orientation
%     clear u s
%     [u, s, v] = svd(squeeze(L2));
%     L_backward(:,is,:) = u(:,:)*s(:,1:2);
% end
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

d=whos; sum([d.bytes])/1000^3
%% calculate MIM
%pca pipeline ('all' 8 pipelines + baseline)
[mic, mim, to_save] = fp_get_mim(A,CS,fqA,nfqA, D,params.ihemi,'all');
fprintf('pipelines calculated')
d=whos; sum([d.bytes])/1000^3

%% performance measures
fprintf('Performance measures... \n')
tic
[PERFORMANCE, BASELINE] = fp_get_performance(gt, mic, mim, params);
toc
d=whos; sum([d.bytes])/1000^3

%% save performance and baseline
fprintf('Saving... \n')
tic
outname = sprintf('%smim_%s.mat',DIROUT,params.logname);
save(outname,'-v7.3')
toc

