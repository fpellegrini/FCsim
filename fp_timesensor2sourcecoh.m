function [P_all, Pv_all, CSv_all, COH_all] = fp_timesensor2sourcecoh(patientNumber)
%pipeline to get from time-series data to coherence on source level

cd ~/Dropbox/MEG_Project/Data

if ~exist('patientNumber','var')
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors
else
    patientID{1} = patientNumber;
end

for id = 1:numel(patientID)
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    
    X = data.D(:,:,:);
    D_ft = ftraw(data.D);
    N_trials = length(D_ft.trial);
    
    fs = data.D.fsample;
    fres = fs;
    frqs = sfreqs(fres, fs);
%     frq_inds = find(frqs > 13 & frqs < 30);
    nfreq = numel(frqs);
    
    id_meg_chan = 1:125;
    id_meg_chan(data.D.badchannels)=[];
    id_lfp_chan = 126:131;
    id_meg_trials = 1:size(X,3);
    id_lfp_trials = id_meg_trials;
    
    %power
    clear Pxx
    for itrial = 1:N_trials
        [Pxx(:, :, itrial), ~] = pwelch(squeeze(X([id_meg_chan id_lfp_chan],:,itrial))', [], [], 2*fs, fs);
    end
    P = mean(Pxx,3);
    
    %cross spectrum
    CS = fp_tsdata_to_cpsd(X,0,fres,id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials); %125x6
    
    %filter
    A1=inverse.MEG.W;
    for i =1:length(A1)
        A(i,:) = real(A1{i});
    end
    
    nvox = size(A,1);
    nlfp = numel(id_lfp_chan);
    %project cross spectrum and power to voxel space
    for ifq = 1:nfreq
        CSv(ifq,:,:) = A * CS(:,:,ifq);
        Pv(ifq,:,:) = cat(2,A,eye(nvox,nlfp)) * P(ifq,:)';
    end
    
    %coherence
    coh = CSv;
    for ijack = 1
        for ifreq = 1:nfreq
            coh(ifreq, :, :, ijack) = squeeze(CSv(ifreq, :, :, ijack)) ...
                ./ sqrt(Pv(ifreq,:)'*P(ifreq,end-nlfp+1:end));
        end
    end
    
    P_all{id} = P;
    Pv_all{id} = Pv;
    CSv_all{id} = CSv;
    COH_all{id} = coh;  
    
    clearvars -except patientID id P_all Pv_all CSv_all COH_all
end

