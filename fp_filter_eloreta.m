function A = fp_filter_eloreta(CS,L)

ns = size(L,2);

for is=1:ns %iterate across nodes 
    Lloc=L(:,is,:);
    filter = squeeze(mkfilt_eloreta_v2(Lloc));

%     for ic = 1:size(Lloc,1)
%         [u,~,~] = svd(squeeze(Lloc(ic,1,:)));
%     end 
        
    %select best orientation 
    csd = filter'*real(CS)*filter;
    [u,~,~] = svd(csd);
    
    A(:,is) = filter*u(:,1);

end