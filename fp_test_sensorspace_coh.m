function fp_test_sensorspace_coh(patientNumber)

cd ~/Dropbox/MEG_Project/Data
DIROUT =  '~/Dropbox/MEG_Project/Data/';

if ~exist('patientNumber','var')
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors
else
    patientID{1} = patientNumber;
end

meg_inds = 1:125;
lfp_inds = 126:131;
N_meg = length(meg_inds);

nlags = 20;
cond = 0;
nboot = 1;
N_it = 100;

for id = 1:numel(patientID)
    
    fileName = sprintf('redPLFP%s_off', patientID{id});
    D = spm_eeg_load(fileName);
    D_ft = ftraw(D);
    N_trials = length(D_ft.trial);
    [~, N_samples] = size(D_ft.trial{1});
    
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frq_inds = find(frqs > 13 & frqs < 30);
    
    data1 = [D_ft.trial{:}];
    
    data1 = reshape(data1([meg_inds lfp_inds], :), N_meg+numel(lfp_inds), N_samples, N_trials);
    
    
    

    for iit = 1: N_it
        
        clear id_lfp_trials
        if iit ==1
            id_lfp_trials = 1: N_trials;
        else
            id_lfp_trials = randperm(N_trials);
        end
               
        data = cat(1,data1(meg_inds,:,:),data1(lfp_inds,:,id_lfp_trials));
        
        conn = data2spwctrgc(data, fres, nlags, cond, nboot, [], {'CS'});
        COH = cs2coh(conn.CS);
        
        r_max_abs(iit,:) = max(mean(mean(abs(COH(frq_inds, meg_inds, 1:3)),1),3));
        r_max_im(iit,:) = max(mean(mean(abs(imag(COH(frq_inds, meg_inds, 1:3))),1),3));
        
        l_max_abs(iit,:) = max(mean(mean(abs(COH(frq_inds, meg_inds, 4:6)),1),3));
        l_max_im(iit,:) = max(mean(mean(abs(imag(COH(frq_inds, meg_inds, 4:6))),1),3));
       
        clear data conn COH 
    end
    
    p_r_abs = sum(r_max_abs(2:end)>r_max_abs(1))/numel(r_max_abs);
    p_l_abs = sum(l_max_abs(2:end)>l_max_abs(1))/numel(l_max_abs);
    p_r_im = sum(r_max_im(2:end)>r_max_im(1))/numel(r_max_im);
    p_l_im = sum(l_max_im(2:end)>l_max_im(1))/numel(l_max_im);
    
    outname = sprintf('%sperm_sensor_Patient%s',DIROUT,patientID{id});
    save(outname,'p_r_abs','p_l_abs','p_r_im','p_l_im','-v7.3')
    clear r_max_abs r_max_im l_max_abs l_max_im p_r_abs p_l_abs p_r_im p_l_im 
    
    
end



