patientNumber = '10';

load(sprintf('Coherences_Patient%s.mat',patientNumber))
[c, voxID] = fp_find_commonvox;

vals = squeeze(abs(coh(1,13,:,1)));
vals1 = vals(voxID{1});


[x y z outtype] = mni2orFROMxyz(c(:,1), c(:,2), c(:,3),[],'mni');



cube = zeros([91 109 91]);

for i =1:numel(x)
    
    cube(x(i),y(i),z(i)) = vals1(i);
end 

cube1 = cube.*100000;
cube2 = round(smooth3(cube1,'box',[11 11 11]));


V= wjn_read_nii('./mri/rPLFP04.nii');

V.fname = 'patient04.nii';

V.img = cube2;
spm_write_vol(V,cube2);