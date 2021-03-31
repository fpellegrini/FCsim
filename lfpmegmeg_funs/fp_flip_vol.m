function [flipped_pos, flip_id] = fp_flip_vol(mni_pos)
%flips positions along sagittal axis and also returns new index

flip_id = nan(size(mni_pos,1),1);
mni_pos= round(mni_pos);

for i = 1:size(mni_pos,1)
    clear a
    a = find(mni_pos(:,1)==-mni_pos(i,1)& mni_pos(:,2)==mni_pos(i,2)...
            & mni_pos(:,3)==mni_pos(i,3));
    if ~isempty(a)
        flip_id(i)= a;
    end
end 

flipped_pos = mni_pos;
flipped_pos(:,1) = -flipped_pos(:,1);