function mni_pos = fp_getMNIpos(patientNumber)
%transforms mni-aligned posisitions to mni positions (even grid with origin
%at zero
% keyboard
load(sprintf('BF_Patient%s.mat',patientNumber))

ori_pos = sources.grid.pos;
trans = data.transforms.toMNI;
mni_pos = cat(2,ori_pos,ones(size(ori_pos,1),1))*trans';
mni_pos = round(mni_pos(:,1:3));