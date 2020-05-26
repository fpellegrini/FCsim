function [nroi,nvox,ind_cortex,ind_roi_cortex, sub_ind_cortex, roi2vox, leadfield] = fp_get_Desikan

load('./processed_bs/bs_results.mat')
% number of ROIs in the Desikan-Kiliany Atlas
nroi = length(cortex.Atlas(3).Scouts);
%roi inds
ind_cortex = [];
sub_ind_cortex = [];
ind_roi = {};
sub_ind_roi = {};
for iROI = 1:nroi
    ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices;
    ind_cortex = cat(1, ind_cortex, ind_roi{iROI});
    [~, ind_roi_cortex{iROI}, ~] = intersect(ind_cortex, ind_roi{iROI}); %index of roi voxels in ind_cortex
    sub_ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices(1);
    sub_ind_cortex = cat(1,sub_ind_cortex, sub_ind_roi{iROI});
    [~,sub_ind_roi_cortex{iROI},~] =  intersect(sub_ind_cortex, sub_ind_roi{iROI});%only one voxel per region
    
end

nvox = length(ind_cortex);
roi2vox = zeros(nvox,1); 
for iroi = 1:nroi 
    roi2vox(ind_cortex(ind_roi_cortex{iroi})) = iroi;
end
roi2vox(roi2vox==0)=[];

