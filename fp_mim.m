function [mic,mim]= fp_mim(Cohroi)

[~,nroi,nfreq] = size(Cohroi,3);
ic=1;
    for iroi = 1:nroi
        
        
        jc=1;
        for jroi = 1:nroi
            
            for ifq = 1:nfreq
                cs_red=[];
                cs_red{1} = Cohroi(ic:ic+npcs-1,ic:ic+npcs-1,ifq);
                cs_red{2} = Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,ifq);
                cs_red{3} = Cohroi(jc:jc+npcs-1,jc:jc+npcs-1,ifq);
                
                caainv=inv(real(cs_red{1})+regu*eye(npcs)*mean(diag(real(cs_red{1}))));
                cab=imag(cs_red{2});
                cbbinv=inv(real(cs_red{3})+regu*eye(npcs)*mean(diag(real(cs_red{3}))));
                X=cab*cbbinv*cab';
                % MIM Ewald Eq. 14
                mim(iroi,jroi,ifq)=(trace(caainv*X));
                caainvsqrt=sqrtm(caainv);
                Y=caainvsqrt*X*caainvsqrt; %Eq. 23
                [~,s,~]=svd(Y);
                % MIC
                mic(iroi,jroi,ifq)=sqrt(s(1,1));
            end
            jc = jc+npcs;
        end
        ic=ic+npcs;
    end