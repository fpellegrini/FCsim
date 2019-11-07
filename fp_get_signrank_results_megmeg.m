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

            [ps(ifreq,iroi,jroi), hs(ifreq,iroi,jroi),stats] = signrank(dbCoh(:,iroi,jroi,ifreq),...
                0,'tail','right','alpha',alpha);
            testval(ifreq,iroi,jroi) = stats.signedrank;
        end
    end
end