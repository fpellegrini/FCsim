function [label,code,id]=fp_get_mni_anatomy(coord)

coord=round(coord);
load('ROI_MNI_V5_List.mat')
t1=wjn_read_nii('ROI_MNI_V5.nii');

[d,xyz]=spm_read_vols(t1);

a = find(xyz(1,:) == coord(1) & xyz(2,:) == coord(2) & xyz(3,:) == coord(3));
id = d(ind2sub(t1.dim,a));
label = [];code =[]; 
if isempty(a)
    id = 0;
% else
%     for i = 1:length(ROI)
%         if ROI(i).ID == id
%             label = ROI(i).Nom_L;
%             code = ROI(i).Nom_C;
%         end
%     end      
end
