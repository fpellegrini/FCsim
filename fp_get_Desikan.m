function D = fp_get_Desikan(iReg)

%load leadfield and cortex structures 
load('processed_bs_wzb_90_2000/bs_results.mat')

%number of ROIs in the Desikan-Kiliany Atlas
nroi = length(cortex.Atlas(3).Scouts);

% all voxel inds
ind_cortex = []; %all voxel indices 
ind_roi = {}; %all voxels grouped by region 
ind_roi_cortex = {}; %index of roi voxels in ind_cortex, grouped by region 

%active voxel inds 
sub_ind_cortex = []; %randomly selected active voxels of each region
sub_ind_roi = {}; %randomly selected active voxels of each region, grouped by region
sub_ind_roi_region = {}; %index of active voxel within region 

for iROI = 1:nroi
    
    ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices;
    ind_cortex = cat(1, ind_cortex, ind_roi{iROI});
    [~, ind_roi_cortex{iROI}, ~] = intersect(ind_cortex, ind_roi{iROI}); 
    sub_ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices(...
        randperm(numel(cortex.Atlas(3).Scouts(iROI).Vertices),iReg));    
    sub_ind_cortex = cat(1,sub_ind_cortex, sub_ind_roi{iROI});
    for ii = 1:iReg
        sub_ind_roi_region{iROI}(ii) = find(ind_roi{iROI}==sub_ind_roi{iROI}(ii));
    end
    
end

%maps roi indeices to voxels 
nvox = length(ind_cortex);
roi2vox = zeros(nvox,1); 
for iroi = 1:nroi 
    roi2vox(ind_cortex(ind_roi_cortex{iroi})) = iroi;
end
roi2vox(roi2vox==0)=[];

%output variable
D.nroi = nroi;
D.nvox = nvox;
D.ind_cortex = ind_cortex;
D.ind_roi_cortex = ind_roi_cortex;
D.sub_ind_cortex = sub_ind_cortex;
D.sub_ind_roi = sub_ind_roi;
D.sub_ind_roi_region = sub_ind_roi_region;
D.roi2vox = roi2vox;
D.leadfield = leadfield;
D.normals = cortex.VertNormals;