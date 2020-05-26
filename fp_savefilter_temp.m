function fp_savefilter_temp(DIROUT)
%pipeline to get from time-series data to coherence on source level

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

for id = [1 2 8]
    
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    
    D=spm_eeg_load(sprintf('redPLFP%s_off',patientID{id}));
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
    id_meg_trials = 1:n_trials;
    id_lfp_trials = 1:n_trials;
    
    %scaling
    X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^6;
    
    %frequency parameters
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frqs(frqs>90) = [];
    nfreq = numel(frqs);
    
    %Now calculate power, cross spectrum and filters. In P, CS, L and A, bad
    %channels are already sorted out.
    
    %cross spectrum
    CS = fp_tsdata_to_cpsd(X,fres,'WELCH',[id_meg_chan,id_lfp_chan], [id_meg_chan,id_lfp_chan], id_meg_trials, id_lfp_trials);
    
    %leadfield
    clear L
    L = fp_get_lf(inverse);
    ns = size(L,2);
    
    %filter
    A=zeros(nmeg,ns,nfreq);
    
    for ifrq = 1:nfreq
        A(:,:,ifrq) = fp_filter(squeeze(CS(1:end-nlfp,1:end-nlfp,ifrq)), L);
    end
    
    CS(:,:,nfreq+1:end)=[];
    
    outname = sprintf('%sFilter_Patient%s_wo_reg',DIROUT, patientID{id});
    save(outname,'A','CS','-v7.3')
    clearvars -except DIROUT DIRLOG logname patientID id
    
    
end

