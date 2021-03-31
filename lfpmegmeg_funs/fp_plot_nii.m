function Y = fp_plot_nii(fileName, dimOrder)
%plots NifTi in 3D space 
%input: fileName of NifTi file, dimOrder: scrolling through third dimension

if nargin < 2
    dimOrder=[1 2 3];
end 

V=spm_vol(sprintf('%s.nii',fileName));
Y=spm_read_vols(V);
Y1 = permute(Y,dimOrder);
imshow3D(Y)
