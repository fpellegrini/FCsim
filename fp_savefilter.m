function fp_savefilter(patientNumber)
%pipeline to get from time-series data to coherence on source level

cd ~/Dropbox/MEG_Project/Data

DIROUT = '~/Dropbox/MEG_Project/Data/';
if ~exist(DIROUT); mkdir(DIROUT); end

if ~exist('patientNumber','var')
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

for id = 1:numel(patientID)
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    
    X = data.D(:,:,:);
    D_ft = ftraw(data.D);
    n_trials = length(D_ft.trial);
    id_meg_trials = 1:n_trials;
    id_lfp_trials = 1:n_trials;
    
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
    
    
   
    %Now calculate power, cross spectrum and filters. In P, CS, L and A, bad
    %channels are already sorted out. 

    %cross spectrum
    %%
    
    id_meg_chan = 1:5;
    CS = fp_tsdata_to_cpsd(X,1,fres,id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials);
    
    %leadfield
    L1 = inverse.MEG.L;
    ns = numel(L1);
    for is=1:ns
        L(:,is,:)= L1{is};
    end     
    
    %filter
    A=zeros(nmeg,ns,nfreq);
    
    for ifrq = 1:nfreq
        currentCS = squeeze(CS(1:end-nlfp,1:end-nlfp,ifrq));
        A(:,:,ifrq) = fp_filter(currentCS, L);
        clear currentCS
    end 
    
    CS(:,:,end-nfreq+1:end)=[];
    
    outname = sprintf('%sFilter_Patient%s',DIROUT, patientID{id});
    save(outname,'A','CS','-v7.3')
    clearvars -except DIROUT patientID id
    
end

