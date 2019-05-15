function A = fp_filter_eloreta(CS,L)

ns = size(L,2);

for is = 1: ns 
    %remove radial orientation 
    clear u v s
    [u v s] = svd(squeeze(L(:,is,:)),'econ');
    L2(:,is,:) = u(:,1:2)*s(1:2,1:2);
end 

filter = squeeze(mkfilt_eloreta_v2(L2));

for is = 1:ns
    
    %select best orientation 
    csd = squeeze(filter(:,is,:))'*real(CS)*squeeze(filter(:,is,:));
    [u,~,~] = svd(csd);
    
    A(:,is) = squeeze(filter(:,is,:))*u(:,1);

end