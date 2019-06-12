function fp_savefilter_eloreta(patientNumber,DIROUT,DIRLOG)
%pipeline to get from time-series data to coherence on source level

fp_addpath

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
        
        load(sprintf('Filter_Patient%s.mat',patientID{id}))
        clear A
        
        load(sprintf('BF_Patient%s.mat',patientID{id}));
        
        L1 = inverse.MEG.L;
        ns = numel(L1);
        for is=1:ns
            L(:,is,:)= L1{is};
        end
        L=L.*(10^(-log10(range(L(:))))); %scale leadfield to avoid numerical issues
        
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
        
        %filter
        A=zeros(nmeg,ns);
        
        for ifrq = 1:nfreq
            currentCS = squeeze(CS(1:end-nlfp,1:end-nlfp,ifrq));
            A(:,:,ifrq) = fp_filter_eloreta(currentCS, L);
            clear currentCS
        end
        
        
        outname = sprintf('%sFilter_Patient%s_e',DIROUT, patientID{id});
        save(outname,'A','CS','-v7.3')
        clearvars -except DIROUT patientID id DIRLOG logname
        
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
    end
    
end

