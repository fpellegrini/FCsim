function fp_megmeg_pipeline

fp_addpath

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

nit= 1000;

for id = 1:numel(patientID)
    
    load(sprintf('Filter_Patient%s.mat',patientID{id}));%Filter and whole CS
    clear CS
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    
    X = D(:,:,:);
    D_ft = ftraw(D);
    n_trials = length(D_ft.trial);
    
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frqs(frqs>90) = [];
    nfreq = numel(frqs);
    
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    id_trials_1 = 1:n_trials;
    
    for iit = 1:nit
        
        if iit ==1
            shuffle = 0;
            id_trials_2 = id_trials_1;
        else
            shuffle =1;
            rng('shuffle')
            id_trials_2 = randperm(n_trials);
        end
        
        
        CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
        
        
        %project cross spectrum to voxel space
        for ifq = 1:nfreq
            CSv(ifq,:,:) = A(:,:,ifq)' * CS(:,:,ifq);
        end
        
        %get voxel power
        for ifq=1:nfreq
            pv(:,ifq) = diag(squeeze(CSv(ifq,:,:)));
        end
        
        %coherence
        coh = CSv;
        for ifreq = 1:nfreq
            coh(ifreq, :, :) = squeeze(CSv(ifreq, :, :)) ...
                ./ sqrt(pv(:,ifreq)*pv(:,ifreq)');
        end
        
        if shuffle == 0
            outname = sprintf('%true_megmeg_coh_Patient%s',DIROUT, patientID{id});
            save(outname,'coh','-v7.3')
        else 
            COH(iit-1,:,:,:) = abs(imag(coh));
        end
        clearvars -except patientID id nit COH A X n_trials nfreq fres n_trials id_meg_chan id_trials_1            
    end
    
    [commonvox_pos, voxID] = fp_find_commonvox;
    mni_pos = fp_getMNIpos(patientID{id});
    conn = fp_find_neighbours(patientID{id}); %think about how to do that
    match_conn = conn(voxID{id},voxID{id});
    ns = size(match_conn,1);
    
    %get flip id and symmetric head
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    match_pos = sym_pos(voxID{id},:);
    [~,flip_id] = fp_flip_vol(match_pos);
    
    freq_conn = zeros(nfreq,nfreq);
    for ifreq = 1:nfreq-1
        freq_conn(ifreq,ifreq+1)=1;
        freq_conn(ifreq+1,ifreq)=1;
    end
    conn_s = sparse(match_conn); %nvox x nvox
    freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
    kron_conn = kron(conn_s,freq_conn_s); %say nvox*nfreq = nkron
    
    clearvars -except patientID id nit 
end

