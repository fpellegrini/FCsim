[pos, voxID] = fp_find_commonvox;

for ii = 1: size(pos,1)
    [label{ii},code{ii},ind(ii)]=fp_get_mni_anatomy(pos(ii,:));
end 

