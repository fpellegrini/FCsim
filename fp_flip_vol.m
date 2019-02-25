function flipped_gridpos = fp_flip_vol(ori_gridpos)

flipped_gridpos = ori_gridpos;
flipped_gridpos(:,1) = -flipped_gridpos(:,1); 