function fp_data2nii(data,pos,outname)

% patientNumber = '10';

% load(sprintf('Coherences_Patient%s.mat',patientNumber));
% vals = squeeze(abs(coh(1,13,:,1)));
% [c, voxID] = fp_find_commonvox;
% data = vals(voxID{1});
% pos = c;

[x y z outtype] = mni2orFROMxyz(pos(:,1), pos(:,2), pos(:,3),[],'mni');
x = round(x);
y=round(y);
z=round(z);

cube = nan([91 109 91]);

for i =1:numel(x)    
    cube(x(i),y(i),z(i)) = data(i);
end 

%  [X,Y,Z] = meshgrid(x1, y1, z1);
%  
% for i=1:numel(x)
%     q(c(i,1)/10+8,c(i,2)/10+12,c(i,3)/10+12)= data(i);
% end 

[X Y Z] = meshgrid(x1,y1,z1);
cube = permute(cube,[2 1 3]);

[x2 y2 z2] = find(isnan(cube));


Vq = interp3(X,Y,Z,cube,x2,y2,z2,'linear');

cube1 = cube.*(11^3);
cube2 = round(smooth3(cube1,'box',[9 9 9]));

V= wjn_read_nii('./mri/rPLFP04.nii');

% V.fname = 'patient04.nii';
V.fname = outname;

V.img = cube2;
spm_write_vol(V,cube2);