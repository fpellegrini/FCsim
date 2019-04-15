

[c, voxID] = fp_find_commonvox;

[x] = min(c);
x=x-1;

y = max(c);

nii_siz = y - x;

Vol = nan(nii_siz);


for is = 1: size(c,1)
    Vol(c(is,1)-x(1),c(is,2)-x(1), c(is,3)-x(3)) = mask(is);
end

nii = make_nii(Vol,[60 100,70]);

niftiwrite(Vol,'nifti3')