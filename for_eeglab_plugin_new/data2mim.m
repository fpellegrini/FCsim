function [MIM,MIC,DIFFGC] = data2mim(signal_sensor,L,ind_roi_cortex)
%Input: 
% signal_sensor in the shape of n_sensor x l_epoch x n_trials
% leadfield L in the shape of n_sensors x n_voxels x 3 dimensions
% ind_roi_cortex which is a 1xnrois cell array that provides the indices of
% the voxels per region 
%Output:
% MIM in the shape of nroi x nroi x nfreq 
% MIC in the shape of nroi x nroi x nfreq
% DIFFGC in the shape of nroi x nroi x nfreq (only upper triangle filled)

% number of sensors, length of epochs, number of trials 
[n_sensors, l_epoch, n_trials] = size(signal_sensor); 
%number of dipole dimensions 
ndim = size(L,3); 
% number of regions of interest 
nroi = numel(ind_roi_cortex);
%frequency resolution 
fres = l_epoch/2; 

%in this configuration npcs (number of principal components) 
%is the same for all regions (and inds too) and is selected according to 
%the minimum number of voxels in a region
for aroi = 1:nroi 
    nvoxroi(aroi) = numel(ind_roi_cortex{aroi});
end
npcs = min(nvoxroi)*ndim;
inds = fp_npcs2inds([npcs npcs]);

%% calculate lcmv for source projection 

%covariance matrix
cCS = cov(signal_sensor(:,:)');
%regularization parameter 
reg = 0.05*trace(cCS)/length(cCS);
%regularized covariance 
Cr = cCS + reg*eye(size(cCS,1));

%calculate lcmv source projection filter A from leadfield L and regularized
%covariance matrix Cr
[~, A] = lcmv(Cr, L, struct('alpha', 0, 'onedim', 0));

%permute A to the shape of n_sensors x ndim x n_voxels
A = permute(A,[1, 3, 2]);

%% PCA

%loop over regions 
for aroi = 1:nroi
    
    clear A_ signal_source
    
    %A_ is the lcmv filter at aroi
    A_ = A(:, :,ind_roi_cortex{aroi},:);
    %number of voxels at the current roi
    nvoxroi(aroi) = size(A_,3); 
    
    A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
    
    %project sensor signal to voxels at the current roi (aroi)
    signal_source = A2{aroi}' * signal_sensor(:,:);
    
    %do PCA 
    clear signal_roi_ S_
    [signal_roi_,~,~] = svd(double(signal_source)','econ');
    
    %bring signal_roi to the shape of npcs x l_epoch x n_trials
    signal_roi{aroi} = reshape(signal_roi_(:,1:npcs)',[],l_epoch,n_trials);
    
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
        [trgc,~,mim,mic, ~, inds] = data2sctrgcmim(data, fres, [], 0, 0, [], inds);
        
        %so far, mim{1} is equal to mim{2} 
        MIM([oroi, uroi],[oroi, uroi],:)=mim{1};
        MIC([oroi, uroi],[oroi, uroi],:)=mic{1};
        DIFFGC(oroi,uroi,:) = squeeze(trgc(:,1) - trgc(:,2));
        
    end
end
toc
