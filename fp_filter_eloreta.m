function A = fp_filter_eloreta(CS,L)
filter = squeeze(mkfilt_eloreta_v2(L));

for is = 1:ns
    
    %select best orientation 
    csd = squeeze(filter(:,is,:))'*real(CS)*squeeze(filter(:,is,:));
    [u,~,~] = svd(csd);
    
    A(:,is) = squeeze(filter(:,is,:))*u(:,1);

end