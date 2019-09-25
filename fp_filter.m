function A = fp_filter(CS, L)
% keyboard
ns = size(L,2);
lambda = mean(diag(real(CS)))/100;

CSinv=pinv(real(CS)); %+lambda * eye(size(CS))

for is=1:ns %iterate across nodes 
    Lloc=squeeze(L(:,is,:));
    filter = pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv; %create filter

    %select best orientation 
    csd = filter*real(CS)*filter';
     
    [u,~,~] = svd(csd);
    LF = Lloc*u(:,1);

    %recompute filter in best orientation 
    A(:,is)=pinv((LF'*CSinv*LF))*LF'*CSinv;

end