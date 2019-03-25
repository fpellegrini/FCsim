function flipped = fp_flip_vol(mni_pos)
%flips positions along sagittal axis 

flipped = mni_pos;
flipped(:,1) = -flipped(:,1); 