function [MIM,MIC,iCOH,DIFFGC,active_rois] = data2mim(signal_sensor,L,ind_roi_cortex)
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
nvoxroi(117)=0;

empty_rois = find(nvoxroi==0);
nvoxroi1 = nvoxroi;
nvoxroi1(empty_rois)=[];

active_rois = 1:nroi;
active_rois(empty_rois)=[];

%% indeces 

% npcs = repmat(min(nvoxroi1)*ndim,1,nroi);
npcs = repmat(min(5,min(nvoxroi1)*ndim),1,nroi);
beg_inds = cumsum([1 npcs(1:end-1)]);
end_inds = cumsum([npcs]);

for iroi = 1:nroi
  PCA_inds{iroi} = beg_inds(iroi):end_inds(iroi);
end

inds = {}; ninds = 0;
for iroi = 1:nroi
    if ismember(iroi,active_rois)
      for jroi = (iroi+1):nroi
          if ismember(jroi,active_rois)
            inds{ninds+1} = {PCA_inds{iroi}, PCA_inds{jroi}};    
        %     inds{ninds+2} = {PCA_inds{jroi}, PCA_inds{iroi}};  
            ninds = ninds + 1;
          end
      end
    end
end


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

signal_roi = []; 
%loop over regions 
for aroi = 1:nroi
    
    if ismember(aroi,active_rois)
        clear A_ signal_source

        %A_ is the lcmv filter at aroi
        A_ = A(:, :,ind_roi_cortex{aroi},:);
        %A_ = A_(:,:,D.sub_ind_roi_region{aroi});
        %number of voxels at the current roi
        nvoxroi(aroi) = size(A_,3);

        A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);

        %project sensor signal to voxels at the current roi (aroi)
        signal_source = A2{aroi}' * signal_sensor(:,:);

        %do PCA
        clear signal_roi_ S_
        [signal_roi_,S_,~] = svd(double(signal_source)','econ');

        signal_roi_ = signal_roi_(:,1:npcs) * S_(1:npcs, 1:npcs);

        %bring signal_roi to the shape of npcs x l_epoch x n_trials
        signal_roi = cat(1,signal_roi,reshape(signal_roi_',[],l_epoch,n_trials));
        %     signal_roi{aroi} = reshape(signal_source,[],l_epoch,n_trials);
    else
        signal_roi = cat(1,signal_roi,zeros(npcs(aroi),l_epoch,n_trials));
    end
    
    
    end



%% calculate MIM/MIC
tic
conn = data2sctrgcmim(signal_roi, fres, 20, 0,0, [], inds, {'MIC', 'MIM','COH'});
toc
%%
nf = size(conn.MIM,1);
MIM = zeros(nroi,nroi,nf);
MIC=zeros(nroi,nroi,nf);
iCOH = zeros(nroi, nroi, nf);

iinds = 0;
for iroi = 1:nroi
    if ismember(iroi,active_rois)
        for jroi = (iroi+1):nroi
            if ismember(jroi,active_rois)
                iinds = iinds + 1; 
                MIM( iroi, jroi,:) = conn.MIM(:, iinds)';
                MIC(iroi, jroi,:) = conn.MIC(:, iinds)';
                iCOH(iroi,jroi,:) = squeeze(mean(mean(abs(imag(conn.COH(:,PCA_inds{iroi},PCA_inds{jroi}))),2),3))';
            end
        end
    end
end

DIFFGC = [];
MIM(empty_rois,:,:)=[];
MIM(:,empty_rois,:)=[];
MIC(empty_rois,:,:)=[];
MIC(:,empty_rois,:)=[];
iCOH(:, empty_rois,:) = []; 
iCOH(empty_rois,:,:)=[]; 


