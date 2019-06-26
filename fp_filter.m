function A = fp_filter(CS, L)

ns = size(L,2);
lambda = mean(diag(real(CS)))/100;

CSinv=pinv(real(CS)+lambda * eye(size(CS)));

for is = 1: ns 
    %remove radial orientation 
    clear u v s
    [u v s] = svd(squeeze(L(:,is,:)),'econ');
    L2(:,is,:) = u(:,1:2)*s(1:2,1:2);
end 

for is=1:ns %iterate across nodes 
    Lloc=squeeze(L2(:,is,:));
    filter = pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv; %create filter

    %select best orientation 
    csd = filter*real(CS)*filter';
     
    [u,~,~] = svd(csd);
    LF = Lloc*u(:,1);

    %recompute filter in best orientation 
    A(:,is)=pinv((LF'*CSinv*LF))*LF'*CSinv;

end