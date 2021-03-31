function [ps,hs,testval,effectdir] = fp_get_signrank_results_gc(DIFFGC,alpha)

[nsubs,nvox,nside,nfreq] = size(DIFFGC);

%sign-rank test
for ifreq = 1:nfreq
    for iside = 1:nside
        for ivox = 1:nvox
            clear ind 
            [ps(ivox,iside,ifreq), hs(ivox,iside,ifreq),stats] = signrank(DIFFGC(:,ivox,iside,ifreq),...
                0,'tail','both','alpha',alpha);
            testval(ivox,iside,ifreq) = stats.signedrank;
            
            [~, ind] = sort(abs(DIFFGC(:,ivox,iside,ifreq)));
            effectdir(ivox,iside,ifreq) = sign(sum(ind.*sign(DIFFGC(ind,ivox,iside,ifreq))));
        end
    end
end