function [COH_all] = fp_timesensor2sourcecoh(patientNumber, shuffle)
%pipeline to get from time-series data to coherence on source level

cd ~/Dropbox/MEG_Project/Data

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

if ~exist('shuffle','var')
    shuffle = 0;
end

for id = 1:numel(patientID)
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    
    X = data.D(:,:,:);
    D_ft = ftraw(data.D);
    n_trials = length(D_ft.trial);
    
    fs = data.D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frqs(frqs>90) = [];
    nfreq = numel(frqs);
    
    id_meg_chan = 1:125;
    id_meg_chan(data.D.badchannels)=[];
    nmeg = numel(id_meg_chan);
    id_lfp_chan = 126:131;
    nlfp = numel(id_lfp_chan);
    id_meg_trials = 1:n_trials;
    
    if shuffle == 1
        id_lfp_trials = randperm(n_trials);
    else
        id_lfp_trials = id_meg_trials;
    end 
    
   
    %Now load/ calculate power, cross spectrum and filters. In P, CS, and A, bad
    %channels are already sorted out.

    %load filter and whole CS
    load(sprintf('Filter_Patient%s.mat',patientID{id}));
    ns = size(A,2);
    
    %throw away all frequencies above 90 Hz (in filters already done) 
    CS(:,:,nfreq+1:end)=[];
    
    %cross spectrum
    if shuffle ==1
        %when trials are shuffled, the CS between meg and lfp must be
        %re-calculated
        cCS = fp_tsdata_to_cpsd(X,0,fres,id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials);
    else 
        cCS = CS(1:(end-nlfp),end-nlfp+1:end,:);
    end
    
    %project cross spectrum to voxel space
    for ifq = 1:nfreq
        CSv(ifq,:,:) = A(:,:,ifq)' * cCS(:,:,ifq);
    end
    
    %get voxel power
    Pv = fp_project_power(CS,A);
   
    %coherence
    coh = CSv;
    for ifreq = 1:nfreq
        clear Plfp
        Plfp = diag(CS(nmeg+1:end,nmeg:end,ifreq)); %power of lfp channels
        coh(ifreq, :, :) = squeeze(CSv(ifreq, :, :)) ...
            ./ sqrt(Pv(:,ifreq)*Plfp');
    end
    
    COH_all{id} = coh;  
    
    if numel(patientID)==1
        clear COH_all
        COH_all = coh;
    end
    
    clearvars -except patientID id P_all Pv_all CSv_all COH_all
end

