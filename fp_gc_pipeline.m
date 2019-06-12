function fp_gc_pipeline

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

[~, voxID] = fp_find_commonvox;
nlags = 20;

for id = 1:numel(patientID)
    
    load(sprintf('Filter_Patient%s.mat',patientID{id}));%load true CS for constructing filters
    clear A
    
    %load data
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    X = D(:,:,:);
    D_ft = ftraw(D);
    n_trials = length(D_ft.trial);
    
    %channel IDs
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    nmeg = numel(id_meg_chan);
    id_lfp_chan = 126:131;
    nlfp = numel(id_lfp_chan);
    
    %frequency parameters
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frqs(frqs>90) = [];
    nfreq = numel(frqs);
    freqs = linspace(0, 1, fres+1);
    maxfreq = fres+1;
    freqs = freqs(1:maxfreq);
    z = exp(-i*pi*freqs);
    
    %construct filters
    
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    L1 = inverse.MEG.L;
    ns_org = numel(L1);
    for is=1:ns_org
        L(:,is,:)= L1{is};
    end
    
    %delete voxels that are not common in all subs
    mni_pos = fp_getMNIpos(patientID{id});
    [~, noEq] = fp_symmetric_vol(mni_pos);
    L(:,noEq,:) = [];
    L = L(:,voxID{id},:);
    ns = numel(voxID{id});
    
    A = nan(3,nmeg,ns,nfreq);
    for ifrq = 1:nfreq
        clear currentCS lambda CSinv
        currentCS = squeeze(CS(1:end-nlfp,1:end-nlfp,ifrq));
        lambda = mean(diag(real(currentCS)))/100;
        CSinv=pinv(real(currentCS)+lambda * eye(size(currentCS)));
        
        for is=1:ns %iterate across nodes
            clear Lloc
            Lloc=squeeze(L(:,is,:));
            A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv); %create filter
        end
    end
    
    % loop over permutations
    for iit = 1:nit
        
        %cross spectrum
        clear CS
        id_trials_1 = 1:n_trials;
        rng('shuffle')
        id_trials_2 = randperm(n_trials);
        CS = fp_tsdata_to_cpsd(X,fres,'MT',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan], id_trials_1, id_trials_2);
        
        %project cross spectrum to voxel space
        cCS = CS(1:(end-nlfp),end-nlfp+1:end,:); %nmeg x nlfp x nfreq
        CSv = zeros(3,ns+nlfp,ns+nlfp,nfreq);
        for ifq = 1:nfreq %%%%%better solution?
            for idir = 1:3
                CSv(idir,1:ns,end-nlfp+1:end,ifq) = squeeze(A(idir,:,:,ifq))' * cCS(:,:,ifq);
                CSv(idir,end-nlfp+1:end,1:ns,ifq)= squeeze(CSv(idir,1:ns,end-nlfp+1:end,ifq))';
                CSv(idir, 1:ns,1:ns,ifq) = squeeze(A(idir,:,:,ifq))' * CS(1:nmeg,1:nmeg,ifq) * squeeze(A(idir,:,:,ifq));
                CSv(idir,end-nlfp+1:end,end-nlfp+1:end,ifq) = CS(end-nlfp+1:end,end-nlfp+1:end,ifq);
            end
        end
        
        G = cpsd_to_autocov(squeeze(CSv(1,:,:,:)), nlags); %%%how to deal here with 3D? Takes long. Possible to select inds first?
        
        inds = {}; ninds = 0;
        for ii = 1:ns % over nvox
            for ij = ns+1:ns+nlfp %over lfps
                inds{ninds+1} = {ii, ij};
                inds{ninds+2} = {ij, ii};
                ninds = ninds + 2;
            end
        end
        
        % (time-reversed) GC just between sender and receiver sets
        
        % loop over sender/receiver combinations to compute time-reversed GC
        for iind = 1:ninds
            if ~isequal(inds{iind}{1}, inds{iind}{2})
                %       disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] -> [' num2str(inds{iind}{2}) ']'])
                
                subset = [inds{iind}{1} inds{iind}{2}];
                nsubsetvars = length(subset);
                subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
                
                % autocovariance to full forward VAR model
                [A1, SIG] = autocov_to_var3(G(subset, subset, :));
                
                % forward VAR model to state space VARMA models
                [eA, eC, eK, eV, eVy] = varma2iss(reshape(A1, nsubsetvars, []), [], SIG, eye(nsubsetvars));
                
                % backward autocovariance to full backward VAR model
                [AR, SIGR] = autocov_to_var3(permute(G(subset, subset, :), [2 1 3]));
                
                % backward VAR to VARMA
                [eAR, eCR, eKR, eVR, eVyR] = varma2iss(reshape(AR, nsubsetvars, []), [], SIGR, eye(nsubsetvars));
                
                % GC and TRGC computation
                GC(:, iind, iit) = iss_SGC(eA, eC, eK, eV, z, subinds{2}, subinds{1});
                GCR = iss_SGC(eAR, eCR, eKR, eVR, z, subinds{2}, subinds{1});
                TRGC(:, iind, iit) = GC(:, iind, iit) - GCR';
            else
                GC(:, iind, iit) = 0;
                TRGC(:, iind, iit) = 0;
            end
        end        
    end 
    
    %%%diff meg-lfp to lfp-meg trgc ? 
    GC_all = GC_all + GC;
    TRGC_all = TRGC_all + TRGC;
    clear GC TRGC
end

