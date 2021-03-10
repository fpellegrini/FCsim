
DIROUT1=[];
params. iReg=1; 
params.iInt = 1;
params.ilag = 2;
params.isnr = 0.5;
params.iss = 0.5;
params.ip=2;
zs=0;

npcs = 1;
inds = fp_npcs2inds([npcs npcs]);

%%

D = fp_get_Desikan(params.iReg);

[sig,brain_noise,sensor_noise, gt,L_save,iroi_seed,iroi_tar,D, fres, n_trials] = fp_generate_mim_signal...
    (params,D,DIROUT1);

%combine noise sources
noise = params.iss*brain_noise + (1-params.iss)*sensor_noise;
noise = noise ./ norm(noise, 'fro');
%combine signal and noise
signal_sensor1 = params.isnr*sig + (1-params.isnr)*noise;
signal_sensor = signal_sensor1 ./ norm(signal_sensor1, 'fro');

%reshape
signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);

n_sensors = size(signal_sensor,1);

%cross spectrum
fprintf('Calculating cross spectrum... \n')
CS = tsdata_to_cpsd_fast(signal_sensor,fres,'WELCH');
CS(:,:,1)=[];
nfreq = size(CS,3);

L=L_save;
L3 = L(:, D.ind_cortex, :);
L_backward = L3;
ni = size(L_backward,3);

L = L_backward;
clear L_backward L_save L3

save('~/Desktop/champ_tests.mat','-v7.3')

%%

%     sigma2_total = mean(diag(cov(signal_sensor(:, :)')));
%     regu = sigma2_total*0.2;
%     sigu = regu*eye(n_sensors);

sigu = cov(sensor_noise');
tic
[~,~,w] = awsm_champ(signal_sensor(:, :), L(:, :) ,...
    sigu, 200, 3, 2, 0);
toc

A = reshape(w',size(L));
A = permute(A,[1, 3, 2]);
A=abs(A);

fqA = ones(1,nfreq);%only one filter for all freqs.
nfqA = 1; 

a = (abs(A)>10^-8);


%% PCA
empty_rois = [];

%loop over regions 
for aroi = 1:D.nroi
    
    clear A_ signal_source
    
    %A_ is the lcmv filter at aroi
    A_ = A(:, :,ind_roi_cortex{aroi},:);
    a_ = a(:,:,ind_roi_cortex{aroi},:);
    %number of voxels at the current roi
    nvoxroi(aroi) = size(A_,3); 
    
    A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
    
    %project sensor signal to voxels at the current roi (aroi)
    signal_source = A2{aroi}' * signal_sensor(:,:);
    
    %sort out the dims and voxels with zero activity 
    s_ = sum(signal_source,2)>10^-8;
    
    %do PCA 
    clear signal_roi_ S_
    [signal_roi_,~,~] = svd(double(signal_source(s_,:))','econ');
    
    %bring signal_roi to the shape of npcs x l_epoch x n_trials
    try
        signal_roi{aroi} = reshape(signal_roi_(:,1:npcs)',[],l_epoch,n_trials);
    catch
        signal_roi{aroi} = zeros(npcs,l_epoch,n_trials);
        empty_rois = [empty_rois aroi];
    end
    
end



%% calculate MIM/MIC

ninds = length(inds); 

%loop over all roi combinations
tic
for oroi = 1:nroi-1
    if ~ismember(oroi,empty_rois)
        for uroi = oroi+1:nroi
            if ~ismember(uroi,empty_rois)
                
                fprintf(['Calculating regions ' num2str(oroi) ' and ' num2str(uroi) '.\n'])
                
                %calculate CS and COH
                clear data CS
                data = cat(1,signal_roi{oroi},signal_roi{uroi});                
                CS = tsdata_to_cpsd_fast(data, fres, 'WELCH', l_epoch);
                
                for ifreq = 1: fres
                    clear pow
                    pow = real(diag(CS(:,:,ifreq)));
                    COH(:,:,ifreq) = CS(:,:,ifreq)./ sqrt(pow*pow');
                end
                
                % loop over sender/receiver combinations to compute MIM
                for iind = 1:ninds
                    if ~isequal(inds{iind}{1}, inds{iind}{2})
                        disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] -> [' num2str(inds{iind}{2}) ']'])
                        
                        %ind configuration
                        subset = [inds{iind}{1} inds{iind}{2}];
                        nsubsetvars = length(subset);
                        subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
                        
                        %MIC and MIM
                        [mic{iind} , mim{iind}] =  roi_mim(COH,subinds);
                    end
                end
                        
                        
                        
%                 clear mim mic trgc -> trgc not working so far! 
%                 [trgc,~,mim,mic, ~, inds] = data2sctrgcmim(data, fres, 20, 0, 0, [], inds);
                
                %so far, mim{1} is equal to mim{2}
                MIM([oroi, uroi],[oroi, uroi],:)=mim{1};
                MIC([oroi, uroi],[oroi, uroi],:)=mic{1};
%                 DIFFGC(oroi,uroi,:) = squeeze(trgc(:,1) - trgc(:,2));
            end
            
        end
    end
end
toc
