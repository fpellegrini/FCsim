function [flip_id, noEq] = fp_get_flip_id(patientNumber,vox_ind)

mni_pos = fp_getMNIpos(patientNumber);
[sym_pos, noEq] = fp_symmetric_vol(mni_pos);
match_pos = sym_pos(vox_ind,:);
[~,flip_id] = fp_flip_vol(match_pos);