
load('mim_pval.mat')

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2
%get rois
clear mni_pos label code roi_id u_roi_id csroi
mni_pos = fp_getMNIpos(patientID{id});
[sym_pos, noEq] = fp_symmetric_vol(mni_pos)
ns = size(mni_pos,1);
for ii = 1: ns
    [label{ii},code{ii},roi_id(ii)]=fp_get_mni_anatomy_new(mni_pos(ii,:));
end
u_roi_id = sort(unique(roi_id));
nroi = numel(u_roi_id)-1;

a = zeros(29,1); 
a(2:29) = sum(true_clu(:,:,4),2);

b = roi_id +1;

c = a(b);

outname = sprintf('mim_cluster.nii');
fp_data2nii(c,mni_pos,[],outname,id)