function fp_savefilter_eloreta_temp_2D(DIROUT)
%pipeline to get from time-series data to coherence on source level

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

for id = 1:numel(patientID)
    fprintf('Working on subject %s',patientID{id})
    
    load(sprintf('Filter_Patient%s.mat',patientID{id}))
    clear A
    
    clear L
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    L = fp_get_lf(inverse);
    ns = size(L,2);
    
    fs = 300;
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
    A = squeeze(mkfilt_eloreta_v2(L));    
    
    outname = sprintf('%sFilter_Patient%s_e2D',DIROUT, patientID{id});
    save(outname,'A','CS','-v7.3')
    clearvars -except DIROUT patientID id DIRLOG logname
    
end

