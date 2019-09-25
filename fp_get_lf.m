function L = fp_get_lf(inverse)
% keyboard
L1 = inverse.MEG.L;
ns = numel(L1);
for is=1:ns
    
    clear L2
    L2 = L1{is}./10^-12;
    
    %remove radial orientation
    clear u s
    [u, s, v] = svd(squeeze(L2));
    L(:,is,:) = u(:,1:2)*s(1:2,1:2);
    
end
