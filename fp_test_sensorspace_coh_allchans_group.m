function fp_test_sensorspace_coh_allchans_group(patientNumber,fband, DIROUT)

fp_addpath

if ~exist(DIROUT); mkdir(DIROUT); end



if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors
else
    patientID{1} = patientNumber;
end

if strcmp(fband,'theta')
    frq_band = [4 8];
elseif strcmp(fband,'alpha')
    frq_band = [7 13];
elseif strcmp(fband,'beta')
    frq_band = [13 30];
elseif strcmp(fband,'gamma_low')
    frq_band = [30 60];
elseif strcmp(fband,'gamma_high')
    frq_band = [60 90];
else 
    warning('Choosing beta frequency band!')
    frq_band = [13 30];
    fband = 'beta';
end

MEG_INDS = 1:125;
LFP_INDS = 126:131;
N_lfp = numel(LFP_INDS);
N_MEG = numel(MEG_INDS);


nlags = 20;
cond = 0;
nboot = 1;
N_it = 1000;
nsubs = numel(patientID);
nfreq = 8;

COH = nan(nsubs,N_it,nfreq,N_MEG, N_lfp);


for id = 1:numel(patientID)
    
    fprintf('Working on subject %s',patientID{id})
    
    fileName = sprintf('redPLFP%s_off', patientID{id});
    D = spm_eeg_load(fileName);
    D_ft = ftraw(D);
    N_trials = length(D_ft.trial);
    [~, N_samples] = size(D_ft.trial{1});
    meg_inds = MEG_INDS;
    meg_inds(D.badchannels) = [];
    N_meg = length(meg_inds);
    
    
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frq_inds = find(frqs > frq_band(1) & frqs < frq_band(2));
    
    data1 = [D_ft.trial{:}];
    
    data1 = reshape(data1([MEG_INDS LFP_INDS], :), N_meg+numel(D.badchannels)+N_lfp, N_samples, N_trials);%bad channels are sorted out later
    
    for iit = 1: N_it
        iit
        clear id_lfp_trials
        if iit==1
            id_lfp_trials = 1: N_trials;
        else
            rng('shuffle')
            id_lfp_trials = randperm(N_trials);
        end
        
        data = cat(1,data1(meg_inds,:,:),data1(LFP_INDS,:,id_lfp_trials));%here bad channels are sorted out
        
        tic
        conn = data2spwctrgc(data, fres, nlags, cond, nboot, [], {'CS'});
        toc
        temp = cs2coh(conn.CS);
        COH(id,iit,:,meg_inds,:) = temp(frq_inds,1:N_meg,end-5:end);        
        clear data conn 
    end        
end



r_abs = squeeze(nanmedian(nanmedian(nansum(abs(COH(:,:,:,:,1:3)),1),3),5));
r_im = squeeze(nanmedian(nanmedian(nansum(abs(imag(COH(:,:,:,:,1:3))),1),3),5));

l_abs = squeeze(nanmedian(nanmedian(nansum(abs(COH(:,:,:,:,4:6)),1),3),5));
l_im = squeeze(nanmedian(nanmedian(nansum(abs(imag(COH(:,:,:,:,4:6))),1),3),5));

for ichan = 1:N_meg
    p_r_abs(ichan) = sum(r_abs(2:end,ichan)>r_abs(1,ichan))/numel(r_abs(:,ichan));
    p_l_abs(ichan) = sum(l_abs(2:end,ichan)>l_abs(1,ichan))/numel(l_abs(:,ichan));
    p_r_im(ichan) = sum(r_im(2:end,ichan)>r_im(1,ichan))/numel(r_im(:,ichan));
    p_l_im(ichan) = sum(l_im(2:end,ichan)>l_im(1,ichan))/numel(l_im(:,ichan));
end

outname = sprintf('%sperm_sensor_%s_group',DIROUT,fband);
save(outname,'p_r_abs','p_l_abs','p_r_im','p_l_im','-v7.3')


