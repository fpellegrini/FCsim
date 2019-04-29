function fp_data2nii(data,pos,outname)

% patientNumber = '10';

% load(sprintf('Coherences_Patient%s.mat',patientNumber));
% vals = squeeze(abs(coh(1,13,:,1)));
% data = vals(voxID{1});
% [c, voxID] = fp_find_commonvox;

[x y z outtype] = mni2orFROMxyz(pos(:,1), pos(:,2), pos(:,3),[],'mni');
x = round(x);
y=round(y);
z=round(z);

cube = zeros([91 109 91]);

for i =1:numel(x)    
    cube(x(i),y(i),z(i)) = data(i);
end 

cube1 = cube.*100000;
cube2 = round(smooth3(cube1,'box',[11 11 11]));

V= wjn_read_nii('./mri/rPLFP04.nii');

% V.fname = 'patient04.nii';
V.fname = outname;

V.img = cube2;
spm_write_vol(V,cube2);