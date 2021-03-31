function [mic,mim]= roi_mim(Cohroi,subinds)

regu = 0.000001; 
[~,~,nfreq] = size(Cohroi);
nroi = length(subinds);

for iroi = 1:nroi
    ipcs = numel(subinds{iroi});

    for jroi = 1:nroi
        jpcs = numel(subinds{jroi});
        
        for ifq = 1:nfreq
            cs_red=[];
            cs_red{1} = Cohroi(subinds{iroi},subinds{iroi},ifq);
            cs_red{2} = Cohroi(subinds{iroi},subinds{jroi},ifq);
            cs_red{3} = Cohroi(subinds{jroi},subinds{jroi},ifq);
            
            caainv=inv(real(cs_red{1})+regu*eye(ipcs)*mean(diag(real(cs_red{1}))));
            cab=imag(cs_red{2});
            cbbinv=inv(real(cs_red{3})+regu*eye(jpcs)*mean(diag(real(cs_red{3}))));
            X=cab*cbbinv*cab';
            % MIM Ewald Eq. 14
            mim(iroi,jroi,ifq)=(trace(caainv*X));
            caainvsqrt=sqrtm(caainv);
            Y=caainvsqrt*X*caainvsqrt; %Eq. 23
            [~,s,~]=svd(Y);
            % MIC
            mic(iroi,jroi,ifq)=sqrt(s(1,1));
        end        
    end
end
