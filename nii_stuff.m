load('p_singlesub_abs.mat')
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
 
 id = 3;
 
 p1 = p{id};
 t1 = TRUE_CLU{id};
 
 mni_pos = fp_getMNIpos(patientID{id});
 
 [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
 
 trans = inv(V.mat);
mni_pos = cat(2,sym_pos,ones(size(sym_pos,1),1))*trans';
mni_pos(:,4)=[];
pos= mni_pos;

 
 
 outname = '~/Desktop/test.nii';
 fp_data2nii(t1,sym_pos,outname)