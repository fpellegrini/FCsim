function fp_savefilter(patientNumber,DIROUT,DIRLOG)
%pipeline to get from time-series data to coherence on source level

fp_addpath

if ~exist(DIROUT); mkdir(DIROUT); end
if ~exist(DIRLOG); mkdir(DIRLOG); end

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

for id = 1:numel(patientID)
    fprintf('Working on subject %s',patientID{id})
    logname = sprintf('%s',patientID{id});
    
    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
        load(sprintf('BF_Patient%s.mat',patientID{id}));
        D=spm_eeg_load(sprintf('redPLFP%s_off',patientID{id}));
        
        X = D(:,:,:);
        X=X./10^(log10(range(X(:)))-2); %scale the data
        D_ft = ftraw(D);
        n_trials = length(D_ft.trial);
        
        fs = D.fsample;
        fres = 75;
        frqs = sfreqs(fres, fs);
        frqs(frqs>90) = [];
        nfreq = numel(frqs);
        
        id_meg_chan = 1:125;
        id_meg_chan(D.badchannels)=[];
        nmeg = numel(id_meg_chan);
        id_lfp_chan = 126:131;
        nlfp = numel(id_lfp_chan);
        id_meg_trials = 1:n_trials;
        id_lfp_trials = 1:n_trials;
        
        
        %Now calculate power, cross spectrum and filters. In P, CS, L and A, bad
        %channels are already sorted out.
        
        %cross spectrum
        
        CS = fp_tsdata_to_cpsd(X,fres,'WELCH',[id_meg_chan,id_lfp_chan], [id_meg_chan,id_lfp_chan], id_meg_trials, id_lfp_trials);
        
        %leadfield
        L1 = inverse.MEG.L;
        ns = numel(L1);
        for is=1:ns
            L(:,is,:)= L1{is};
        end
        L=L.*(10^(-log10(range(L(:))))); %scale leadfield to avoid numerical issues
        
        %filter
        A=zeros(nmeg,ns,nfreq);
        
        for ifrq = 1:nfreq
            currentCS = squeeze(CS(1:end-nlfp,1:end-nlfp,ifrq));
            A(:,:,ifrq) = fp_filter(currentCS, L);
            clear currentCS
        end
        
        CS(:,:,nfreq+1:end)=[];
        
        outname = sprintf('%sFilter_Patient%s',DIROUT, patientID{id});
        save(outname,'A','CS','-v7.3')
        clearvars -except DIROUT DIRLOG logname patientID id
        
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
    end
    
end

