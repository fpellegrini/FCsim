function [ps,hs,testval] = fp_get_signrank_results_megmeg(tCoh,sCoh,alpha)
%tCoh = coherence to be tested, sCoh = shuffled coherences 
%debias
o = squeeze(median(sCoh,2)); 
dbCoh = tCoh-o; 

[~,nroi,nroi,nfreq] = size(dbCoh);

%sign-rank test
for ifreq = 1:nfreq
    for iroi = 1:nroi
        for jroi = 1:nroi

            [ps(iroi,jroi,ifreq), hs(iroi,jroi,ifreq),stats] = signrank(dbCoh(:,iroi,jroi,ifreq),...
                0,'tail','right','alpha',alpha);
            testval(iroi,jroi,ifreq) = stats.signedrank;
        end
    end
end