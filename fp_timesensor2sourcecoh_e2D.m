function [COH_all] = fp_timesensor2sourcecoh_e2D(patientNumber, shuffle)
%pipeline to get from time-series data to coherence on source level

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

if ~exist('shuffle','var')
    shuffle = 0;
end

for id = 1:numel(patientID)
%     load(sprintf('BF_Patient%s.mat',patientID{id}));
    load(sprintf('Filter_Patient%s_e2D.mat',patientID{id}));%Filter and whole CS
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    
    X = D(:,:,:);
    D_ft = ftraw(D);
    n_trials = length(D_ft.trial);
    
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frqs(frqs>90) = [];
    nfreq = numel(frqs);
    ns = size(A,2);
    
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    nmeg = numel(id_meg_chan);
    id_lfp_chan = 126:131;
    nlfp = numel(id_lfp_chan);
    id_meg_trials = 1:n_trials;
    ndim = size(A,3);
    
    %scaling
    X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./(10^6);
    
    if shuffle == 1
        rng('shuffle')
        id_lfp_trials = randperm(n_trials);
    else
        id_lfp_trials = id_meg_trials;
    end 
    
    %throw away all frequencies above 90 Hz (in filters already done) 
    CS(:,:,nfreq+1:end)=[];
    
    %cross spectrum
    if shuffle ==1
        %when trials are shuffled, the CS between meg and lfp must be
        %re-calculated
        tic
        cCS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials);
        toc
    else 
        cCS = CS(1:(end-nlfp),end-nlfp+1:end,:);
    end
    
    %project cross spectrum to voxel space
    for idim =1:ndim
        for ifq = 1:nfreq
            CSv(ifq,idim,:,:) = squeeze(A(:,:,idim))' * cCS(:,:,ifq);
        end
         
        %get voxel power
        pv(idim,:,:) = fp_project_power_2D(CS(1:nmeg,1:nmeg,:),squeeze(A(:,:,idim)))';
        
    end
   
%     %sum across dimensions 
%     CSv = squeeze(sum(CSv,2));
%     pv = squeeze(sum(pv,1))';
    
    %coherence
    coh = CSv;
    for idim = 1:ndim
        for ifreq = 1:nfreq
            clear plfp
            plfp = diag(squeeze(CS(nmeg+1:end,nmeg+1:end,ifreq))); %power of lfp channels
            coh(ifreq, idim,:, :) = squeeze(CSv(ifreq, idim,:, :)) ...
                ./ sqrt(squeeze(pv(idim,ifreq,:))*plfp');
        end
    end
    
    coh = squeeze(sum(coh,2));
    
    
    COH_all{id} = coh;
    
    if numel(patientID)==1
        clear COH_all
        COH_all = coh;
    end
    
    clearvars -except patientID id COH_all
end

