

clear all
load('BF_Patient08.mat')
pos= sources.pos;
d = eucl(pos,pos(2000,:));
xx = exp(-10^-1.5*d);

scatter3(pos(:,1),pos(:,2),pos(:,3),10,xx)
 
outname = 'xx.nii';
fp_data2nii(xx,pos,[],outname)


%% 
id = 3;
[commonvox_pos, voxID] = fp_find_commonvox;

pos1 = commonvox_pos; 
mni_pos = fp_getMNIpos('08');
[sym_pos, noEq] = fp_symmetric_vol(mni_pos);
pos2 = sym_pos(voxID{id},:);


xx1 = xx;
xx1(noEq)=[];
xx1=xx1(voxID{id});
scatter3(commonvox_pos(:,1),commonvox_pos(:,2),commonvox_pos(:,3),10,xx1)