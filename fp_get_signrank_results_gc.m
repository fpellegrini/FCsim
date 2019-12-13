function [ps,hs,testval] = fp_get_signrank_results_gc(DIFFGC,alpha)

[nsubs,nvox,nside,nfreq] = size(DIFFGC);

%sign-rank test
for ifreq = 1:nfreq
    for iside = 1:nside
    for ivox = 1:nvox
        
        [ps(ivox,iside,ifreq), hs(ivox,iside,ifreq),stats] = signrank(DIFFGC(:,ivox,iside,ifreq),...
            0,'tail','right','alpha',alpha);
        testval(ivox,iside,ifreq) = stats.signedrank;
    end
    end
end