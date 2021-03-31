function [gc, trgc, diffgc] = fp_gc(CSroi,npcs,z)

nlags= 20; 

clear cCS CSinv currentCS G
G = cpsd_to_autocov(CSroi, nlags);

% indices of required channel combinations
inds=[]; ninds=0;
c1 = 0;
for ipcs = 1:numel(npcs)
    c2 = c1+npcs(ipcs);
    for jpcs = (ipcs+1):numel(npcs)
        inds{ninds+1} = {[c1+1:c1+npcs(ipcs)] , [c2+1:c2+npcs(jpcs)]};
        inds{ninds+2} = {[c2+1:c2+npcs(jpcs)], [c1+1:c1+npcs(ipcs)]};
        c2 = c2+npcs(jpcs);
        ninds= ninds+2;
    end
    c1 = c1+npcs(ipcs);
end


%%
% loop over sender/receiver combinations to compute time-reversed GC
for iind = 1:ninds
    %       disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] -> [' num2str(inds{iind}{2}) ']'])
    clear subset subsetvars subinds A1 SIG eA eC eK eV AR SIGR eAR eCR eKR eVR GCR
    subset = [inds{iind}{1} inds{iind}{2}];
    nsubsetvars = length(subset);
    subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
    
    % autocovariance to full forward VAR model
    [A1, SIG] = autocov_to_var4(G(subset, subset, :));
    
    % forward VAR model to state space VARMA models
    [eA, eC, eK, eV, ~] = varma2iss(reshape(A1, nsubsetvars, []), [], SIG, eye(nsubsetvars));
    
    % backward autocovariance to full backward VAR model
    [AR, SIGR] = autocov_to_var4(permute(G(subset, subset, :), [2 1 3]));
    
    % backward VAR to VARMA
    [eAR, eCR, eKR, eVR, ~] = varma2iss(reshape(AR, nsubsetvars, []), [], SIGR, eye(nsubsetvars));
    
    % GC and TRGC computation
    gc(iind,:) = iss_SGC(eA, eC, eK, eV, z, subinds{2}, subinds{1});
    gcr = iss_SGC(eAR, eCR, eKR, eVR, z, subinds{2}, subinds{1});
    trgc(iind,:) = gc(iind,:) - gcr;
end

%%
%diff
o=1;
for iind = 1:2:ninds-1
    diffgc(o,:) = trgc(iind+1,:)-trgc(iind,:);
    o=o+1;
end



