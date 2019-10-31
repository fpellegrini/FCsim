function [ps,hs,testval] = fp_get_signrank_results(tCoh,sCoh,alpha)
%tCoh = coherence to be tested, sCoh = shuffled coherences 

%debias
o = squeeze(median(sCoh,2)); 
dbCoh = tCoh-o; 

[~,nfreq,ns] = size(dbCoh);

%sign-rank test
for ifreq = 1:nfreq
    for ivox = 1:ns
%         [hn(ivox,ifreq), pn(ivox,ifreq)] = kstest(dbCoh(:,ifreq,ivox)); %not one is n.d.
%         [ht(ivox,ifreq), pt(ivox,ifreq),stats] = ttest(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',0.001);
%         testval = stats(1);
        [ps(ifreq,ivox), hs(ifreq,ivox),stats] = signrank(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',alpha);
        testval(ifreq,ivox) = stats.signedrank;
    end
end