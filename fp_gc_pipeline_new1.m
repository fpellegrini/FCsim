function fp_gc_pipeline_new1(patientNumber, DIROUT)

% fp_addpath_sabzi

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

if ~exist('DIROUT','var')
    error('Please indicate where the results should be saved.')
end

nit=10; %%%%
[~, voxID] = fp_find_commonvox;
nlags = 20;
ndim = 2;

nsub = numel(patientID);
nvox = 10; 
nfreq= 46; 

GC= nan(nsub,nvox,4,nfreq);
TRGC= nan(nsub,nvox,4,nfreq);
DIFFGC= nan(nsub,nvox,2,nfreq);


%%

for id = 1:numel(patientID)
    fprintf('Working on subject %d. \n',id)
    %load data
    clear X
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    X = D(:,:,:);
    D_ft = ftraw(D);
    n_trials = length(D_ft.trial);
    
    %channel IDs
    clear id_meg_chan id_lfp_chan
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    nmeg = numel(id_meg_chan);
    id_lfp_chan = 126:131;
    nlfp = numel(id_lfp_chan);
    
    %scaling
    load('scaling_factor.mat')
    X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./(10^6);
    
    %frequency parameters
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    maxfreq = 90;
    frqs(frqs>maxfreq) = [];
    nfreq = numel(frqs);
    z = exp(-i*pi*(frqs./(fs/2)))';
%     z = (-i*pi*frqs)';
    
    %construct filters
    
    %leadfield
    clear L
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    L = fp_get_lf(inverse);
    
    %delete voxels that are not common in all subs
    clear mni_pos noEq
    mni_pos = fp_getMNIpos(patientID{id});
    [~, noEq] = fp_symmetric_vol(mni_pos);
    L(:,noEq,:) = [];
    L = L(:,voxID{id},:);
    ns = numel(voxID{id});
    
    % true CS
    load(sprintf('Filter_Patient%s_e.mat',patientID{id}));% 1D-A and true CS
    clear A
    %     clear id_trials_1 id_trials_2 A_ CS
    %     id_trials_1 = 1:n_trials;
    %     id_trials_2 = 1:n_trials;
    %     CS = fp_tsdata_to_cpsd(X,fres,'MT',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan], id_trials_1, id_trials_2);
    
    %construct filter
    A = squeeze(mkfilt_eloreta_v2(L));
    A = permute(A,[1, 3, 2]);
    
    for ivox=1:10
        tic
        A_ = squeeze(A(:,:,ivox));
        
        %% GC
        
        %project cross spectrum to voxel space
        clear CSv
        cCS = CS(1:(end-nlfp),end-nlfp+1:end,:); %nmeg x nlfp x nfreq
        CSv = zeros(ndim+(nlfp),ndim+(nlfp),nfreq);
        for ifq = 1:nfreq
            
            csv = zeros(ndim+(nlfp),ndim+(nlfp));
            csv(1:ndim,1:ndim) = A_' * CS(1:(end-nlfp),1:(end-nlfp),ifq) * A_;
            csv(1:ndim,end-(nlfp)+1:end) = A_' * cCS(:,:,ifq); %meg lfp
            csv(end-nlfp+1:end,1:ndim)= csv(1:ndim,end-nlfp+1:end)'; %lfp meg
            csv(end-nlfp+1:end,end-nlfp+1:end) = CS(end-nlfp+1:end,end-nlfp+1:end,ifq); %lfp lfp
            
            %replace power with real values
            clear n
            n = size(csv,1);
            csv(1:(n+1):end) = real(diag(csv));
            
            CSv(:,:,ifq) = csv;
            clear csv
        end
        
        clear cCS CSinv currentCS G
        G = cpsd_to_autocov(CSv, nlags);
        
        %%
        
        %indices of required channel combinations
        clear inds
        inds = {}; ninds = 0;
        inds{1} = {[1,2], [3,4,5]}; %two dimensions of the voxel x 3 left lfp channels 
        inds{2} = {[3,4,5],[1,2]}; % 3 left lfp channels x two dimensions of the voxel 
        inds{3} = {[1,2], [6,7,8]}; %two dimensions of the voxel x 3 right lfp channels 
        inds{4} = {[6,7,8],[1,2]}; % 3 right lfp channels x two dimensions of the voxel
        ninds=4;
        
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
            GC(id, ivox,iind,:) = iss_SGC(eA, eC, eK, eV, z, subinds{2}, subinds{1});
            GCR = iss_SGC(eAR, eCR, eKR, eVR, z, subinds{2}, subinds{1});
            TRGC(id,ivox,iind,:) = squeeze(GC(id, ivox,iind,:)) - GCR';
        end
        
        %diff meg-lfp to lfp-meg trgc
        o=1;
        for iind = 1:2:ninds-1
            DIFFGC(id,ivox,o,:) = TRGC(id,ivox,iind+1,:)-TRGC(id,ivox,iind,:);
            o=o+1;
        end
        
        toc
    end
end

outname = sprintf('%sDIFFGC',DIROUT);
save(outname,'GC','TRGC','DIFFGC','-v7.3')

