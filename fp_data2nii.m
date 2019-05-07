function fp_data2nii(data,pos,outname)

% patientNumber = '10';

% load(sprintf('Coherences_Patient%s.mat',patientNumber));
% vals = squeeze(abs(coh(1,13,:,1)));
% [c, voxID] = fp_find_commonvox;
% data = vals(voxID{1});
% pos = c;

mask= wjn_read_nii('/Users/franziskapellegrini/Dropbox/Master/Masterarbeit/MasterThesis/wjn_toolbox/spm12/toolbox/FieldMap/brainmask.nii');
mask = mask.img;

cube = nan([91 109 91]);
cube(mask==0) = 0;

[x y z outtype] = mni2orFROMxyz(pos(:,1), pos(:,2), pos(:,3),[],'mni');
x = round(x);
y=round(y);
z=round(z);

for i =1:numel(x)    
    cube(x(i),y(i),z(i)) = data(i);
end 
Vq = inpaintn(cube,500);
Vq = Vq.*(10^3);
% Vq = inpaintn(cube);


V= wjn_read_nii('./mri/rPLFP04.nii');

% V.fname = 'patient04.nii';
V.fname = outname;

V.img = Vq;
spm_write_vol(V,Vq);