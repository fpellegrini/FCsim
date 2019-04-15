function fp_test_sensorspace_coh_allchans(patientNumber,DIROUT)

if ~exist('DIROUT','var')
    DIROUT =  '~/Dropbox/Data_MEG_Project/';
end
if ~exist(DIROUT); mkdir(DIROUT); end

if ~exist('patientNumber','var')
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors
else
    patientID{1} = patientNumber;
end

MEG_INDS = 1:125;
LFP_INDS = 126:131;
N_lfp = numel(LFP_INDS);


nlags = 20;
cond = 0;
nboot = 1;
N_it = 1000;

for id = 1:numel(patientID)
    
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
    frq_inds = find(frqs > 13 & frqs < 30);
    
    data1 = [D_ft.trial{:}];
    
    data1 = reshape(data1([MEG_INDS LFP_INDS], :), N_meg+numel(D.badchannels)+N_lfp, N_samples, N_trials);%bad channels are sorted out later 

    for iit = 1: N_it
        
        clear id_lfp_trials
        if iit ==1
            id_lfp_trials = 1: N_trials;
        else
            id_lfp_trials = randperm(N_trials);
        end
               
        data = cat(1,data1(meg_inds,:,:),data1(LFP_INDS,:,id_lfp_trials));%here bad channels are sorted out 
        
        conn = data2spwctrgc(data, fres, nlags, cond, nboot, [], {'CS'});
        COH = cs2coh(conn.CS);
        
        r_max_abs(iit,:) = mean(mean(abs(COH(frq_inds,1:N_meg, 1:3)),1),3);
        r_max_im(iit,:) = mean(mean(abs(imag(COH(frq_inds,1:N_meg, 1:3))),1),3);

        l_max_abs(iit,:) = mean(mean(abs(COH(frq_inds, 1:N_meg, 4:6)),1),3);
        l_max_im(iit,:) = mean(mean(abs(imag(COH(frq_inds, 1:N_meg, 4:6))),1),3);
       
        clear data conn COH 
    end
    
    for ichan = 1:N_meg
        p_r_abs(ichan) = sum(r_max_abs(2:end,ichan)>r_max_abs(1,ichan))/numel(r_max_abs(:,ichan));
        p_l_abs(ichan) = sum(l_max_abs(2:end,ichan)>l_max_abs(1,ichan))/numel(l_max_abs(:,ichan));
        p_r_im(ichan) = sum(r_max_im(2:end,ichan)>r_max_im(1,ichan))/numel(r_max_im(:,ichan));
        p_l_im(ichan) = sum(l_max_im(2:end,ichan)>l_max_im(1,ichan))/numel(l_max_im(:,ichan));
    end
    
    outname = sprintf('%sperm_sensor_allchans_Patient%s',DIROUT,patientID{id});
    save(outname,'p_r_abs','p_l_abs','p_r_im','p_l_im','-v7.3')
    clear r_max_abs r_max_im l_max_abs l_max_im p_r_abs p_l_abs p_r_im p_l_im 
    
    
end



