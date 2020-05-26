function [GC, TRGC] = fp_shuffledGC1(z,nmeg,n_trials,X,fres,id_meg_chan,id_lfp_chan,nlfp,nfreq,ndim,ns,A_,nlags)

try
    %cross spectrum
    clear CS id_trials_1 id_trials_2
    id_trials_1 = 1:n_trials;
    rng('shuffle')
    id_trials_2 = randperm(n_trials);
    CS = fp_tsdata_to_cpsd(X,fres,'WELCH',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan], id_trials_1, id_trials_2);
    
    %project cross spectrum to voxel space
    clear CSv
    cCS = CS(1:(end-nlfp),end-nlfp+1:end,:); %nmeg x nlfp x nfreq
    CSv = zeros(ndim*ns+nlfp,ndim*ns+nlfp,nfreq);
    for ifq = 1:nfreq
        
        csv = zeros(ns*ndim+nlfp,ns*ndim+nlfp);
        csv(1:ns*ndim,end-nlfp+1:end) = squeeze(A_(:,:,ifq))' * cCS(:,:,ifq); %meg lfp
        csv(end-nlfp+1:end,1:ns*ndim)= csv(1:ns*ndim,end-nlfp+1:end)'; %lfp meg
        csv(end-nlfp+1:end,end-nlfp+1:end) = CS(end-nlfp+1:end,end-nlfp+1:end,ifq); %lfp lfp
        
        csv1 = squeeze(A_(:,:,ifq))' * CS(1:nmeg,1:nmeg,ifq) * squeeze(A_(:,:,ifq)); %megmeg
        pv = fp_project_power(CS(1:nmeg,1:nmeg,ifq),A_(:,:,ifq)); %voxel power
        n1 = size(csv1,1);
        csv1(1:(n1+1):end) = pv;
        csv(1:ns*ndim,1:ns*ndim)=csv1;
        clear csv1 pv n1
        
        %replace power with real values
        clear n
        n = size(csv,1);
        csv(1:(n+1):end) = real(diag(csv));
        CSv(:,:,ifq) = csv; %.*10^4; %re-scale to avoid numerical errors
        clear csv
        
    end
    clear cCS CS CSinv currentCS G inds
    
    G = cpsd_to_autocov(CSv, nlags);
    
    inds = {}; ninds = 0;
    for ii = 1:ndim:ns*ndim % over nvox
        inds{ninds+1} = {[ii, ii+1], [ns*ndim+1:ns*ndim+nlfp]};
        inds{ninds+2} = {[ns*ndim+1:ns*ndim+nlfp], [ii, ii+1]};
        ninds = ninds + 2;
    end
    
    % (time-reversed) GC just between sender and receiver sets
    
    % loop over sender/receiver combinations to compute time-reversed GC
    for iind = 1:ninds
        if ~isequal(inds{iind}{1}, inds{iind}{2})
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
            GC(:, iind) = iss_SGC(eA, eC, eK, eV, z, subinds{2}, subinds{1});
            GCR = iss_SGC(eAR, eCR, eKR, eVR, z, subinds{2}, subinds{1});
            TRGC(:, iind) = GC(:, iind) - GCR';
        else
            GC(:, iind) = 0;
            TRGC(:, iind) = 0;
        end
    end
catch
    warning('Repeating iteration!')
    [GC, TRGC] = fp_shuffledGC(z,nmeg,n_trials,X,fres,id_meg_chan,id_lfp_chan,nlfp,nfreq,ndim,ns,A_,nlags);
end