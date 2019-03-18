function [P_all, Pv_all, CSv_all, COH_all] = fp_timesensor2sourcecoh(patientNumber, shuffle)
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
    fres = fs;
    frqs = sfreqs(fres, fs);
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
    
   
    %Now calculate power, cross spectrum and filters. In P, CS, L and A, bad
    %channels are already sorted out. 
    
    %power
    P = fp_get_pow(X, fres, id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials);
    
    %cross spectrum
    CS = fp_tsdata_to_cpsd(X,0,fres,id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials);
    
    %leadfield
    L1 = inverse.MEG.L;
    ns = numel(L1);
    for is=1:ns
        L(:,is,:)= L1{is};
    end     
    
    %filter
    load(sprintf('Filter_Patient%s.mat',patientID{id}));

%     A1=inverse.MEG.W;
%     for i =1:length(A1)
%         A2(i,:) = real(A1{i});
%     end
    
    %project cross spectrum and power to voxel space
    for ifq = 1:nfreq
        CSv(ifq,:,:) = A(:,:,ifq)' * CS(:,:,ifq);
        Pv(ifq,:,:) = cat(2,A(:,:,ifq)',eye(ns,nlfp)) * P(:,ifq);
    end
    
    %coherence
    coh = CSv;
    for ifreq = 1:nfreq
        coh(ifreq, :, :) = squeeze(CSv(ifreq, :, :)) ...
            ./ sqrt(Pv(ifreq,:)'*P(end-nlfp+1:end,ifreq)');
    end
    
    P_all{id} = P;
    Pv_all{id} = Pv;
    CSv_all{id} = CSv;
    COH_all{id} = coh;  
    
    if numel(patientID)==1
        clear COH_all
        COH_all = coh;
    end
    
    clearvars -except patientID id P_all Pv_all CSv_all COH_all
end

