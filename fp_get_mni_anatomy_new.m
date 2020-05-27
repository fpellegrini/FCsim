function [label,code,id_new,partner_rois]=fp_get_mni_anatomy_new(coord)
%keyboard
coord=round(coord);
load('ROI_MNI_V5_List.mat')
t1=wjn_read_nii('ROI_MNI_V5.nii');
load('ROI_new_j.mat')
nroi = length(ROI_new.label);

[d,xyz]=spm_read_vols(t1);

a = find(xyz(1,:) == coord(1) & xyz(2,:) == coord(2) & xyz(3,:) == coord(3));

if isempty(a)
    id_new = 0;
else
    id_old = d(ind2sub(t1.dim,a));
    
    if id_old == 0
        id_new = 0;
    else
        for ii =1:length(ROI_new.id_old)
            if ismember(id_old,ROI_new.id_old{ii})
                id_new = ROI_new.id_new(ii);
                break
            end
        end
        if ~exist('id_new')
            id_new = 0;
        end
    end
    
end

label=[];
code =[];

%respective rois on the other hemisphere 
partner_rois(1,:)= 1:nroi; 
partner_rois(2,:) = [15:26, nan, nan,1:12, 28,27];