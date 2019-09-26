function fp_data2nii(data,id, nit, outname)
%data should have be of size 1 x nvox 
%id is index of patient in patientID variable, if data is group data id =
%nan.
%nit is number of iterations (default: 500)
%outname example: 'test.nii'
keyboard
if isempty(nit)
    nit=50;
end

[pos, voxID] = fp_find_commonvox;
if ~isnan(id)
    data = data(voxID{id}); %start data
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

x(z<0)=[];
y(z<0)=[];
z(z<0)=[];

%fill data into cube
for i =1:numel(x)    
    cube(x(i),y(i),z(i)) = data(i);
end 

%interpolation of missing values
Vq = cube.*(10^3);
Vq = inpaintn(Vq,nit);

%read in some template 
V= wjn_read_nii('./mri/rPLFP04.nii');
V.fname = outname;
V.img = Vq;

%write nifti 
spm_write_vol(V,Vq);