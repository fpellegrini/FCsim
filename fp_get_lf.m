function L = fp_get_lf(inverse)

L1 = inverse.MEG.L;
ns = numel(L1);
for is=1:ns
    
    L2(:,is,:)= L1{is};
    
    %remove radial orientation
    clear u s
    [u, ~, s] = svd(squeeze(L2(:,is,:)),'econ');
    L(:,is,:) = u(:,1:2)*s(1:2,1:2);
    
end
L = L./10^-12;