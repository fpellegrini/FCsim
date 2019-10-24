function fp_data2nii(data, pos, nit, outname)
%data should have be of size 1 x nvox 
%pos are the positions in mni format 
%nit is number of iterations (default: 500)
%outname example: 'test.nii'
% keyboard

sf = 10^6;%scaling factor 

if isempty(nit)
    nit=50;
end

if sum(data)==0
    error('data is all zero') 
end
cube = nan([91 109 91]); %destination cube

%mask for all data points outside the brain
mask= wjn_read_nii('/Users/franziskapellegrini/Documents/Master/Masterarbeit/MasterThesis/wjn_toolbox/spm12/toolbox/FieldMap/brainmask.nii');
mask = mask.img;
cube(mask==0) = 0;

%from mni to world coordinates
[x, y, z, ~] = mni2orFROMxyz(pos(:,1), pos(:,2), pos(:,3),[],'mni');
x = round(x);
y=round(y);
z=round(z);

data(z<1) = [];
x(z<1)=[];
y(z<1)=[];
z(z<1)=[];


%fill data into cube
for i =1:numel(x)    
    cube(x(i),y(i),z(i)) = data(i);
end 

%scale
Vq = cube.*sf;

% Vq = inpaintn(Vq,nit);

%read in some template 
V= wjn_read_nii('./mri/rPLFP22.nii');
V.fname = outname;
V.img = Vq;

%write nifti 
spm_write_vol(V,Vq);

X = wjn_read_nii(['./', outname,]);
s=20;
wjn_nii_smooth(X.fname,s)