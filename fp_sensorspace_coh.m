function fp_sensorspace_coh(patientNumber)

cd ~/Dropbox/MEG_Project/Data
DIROUT = '~/Dropbox/MEG_Project/Data/figures/sensorspace_coh/';
if ~exist(DIROUT); mkdir(DIROUT); end

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

for id = 1:numel(patientID)
    
    fileName = sprintf('redPLFP%s_off', patientID{id});
    D = spm_eeg_load(fileName);
    D_ft = ftraw(D);
    N_trials = length(D_ft.trial);
    [~, N_samples] = size(D_ft.trial{1});

    fs = D.fsample;
    fres = fs;
    frqs = sfreqs(fres, fs);
    
    data = [D_ft.trial{:}];
        
    data = reshape(data([meg_inds lfp_inds], :), N_meg+numel(lfp_inds), N_samples, N_trials);

    conn = data2spwctrgc(data, fres, nlags, cond, nboot, [], {'CS'});        
    COH = cs2coh(conn.CS);
    absCOH = abs(COH(:, meg_inds, lfp_inds));
    
    loc = mk_sensors_plane(D_ft.grad.chanpos(:, [2 1 3]));
    close all
    frq_inds = find(frqs > 15 & frqs < 23);
    
%     figone(17,38)
%     for i = 1:6
%         subplot(2,3,i)
%         showfield_general(mean(mean(absCOH(frq_inds, :,i),1),3)', loc); %mean across frequencies of interest and 6 lfp chans
%         caxis([0 0.15])
%         % colormap jet
%         colorbar    
%     end
%     outname = sprintf('%ssix_lfps_Patient%s.png',DIROUT,patientID{id});
%     print(outname,'-dpng');
%     close all

    %right lfp channels
    
    r_COH = mean(absCOH(frq_inds,:,1:3),3); %15 x 125
    temp = max(max(r_COH));
    [r_max_f r_max_chan] = find(r_COH==temp);
    clear temp
    
    figone(30,30)
    pars.scale = [0 0.15];
    %plot right channels
    %mean across freqs between 15 and 22
    subplot(2,2,1)
    showfield_general(mean(r_COH,1)', loc, pars); 
    caxis([0 0.15])
    text(-0.5, 0.7,'mean coherence across frequencies (15-22 Hz)')
    
    %only for the strongest frequency
    subplot(2,2,2)
    showfield_general(r_COH(r_max_f, :)', loc, pars);
    text(-0.4, 0.7,'coherency at strongest frequency')
    
    %coherence spectrum at strongest channel
    subplot(2,2,3)
    plot(frqs(frq_inds), r_COH(:,r_max_chan))
    xlabel('Freq (Hz)')
    ylabel('Coherence')
    title('Coherence spectrum at strongest channel') 
    
    clear outname
    outname = sprintf('%scoh_analysis_Patient%s_right.png',DIROUT,patientID{id});
    print(outname,'-dpng');
    close all
    
    
    %left lfp channels
    
    l_COH = mean(absCOH(frq_inds,:,4:6),3); %15 x 125
    temp = max(max(l_COH));
    [l_max_f l_max_chan] = find(l_COH==temp);
    clear temp
    
    figone(30,30)
    pars.scale = [0 0.15];
    %plot right channels
    %mean across freqs between 15 and 22
    subplot(2,2,1)
    showfield_general(mean(l_COH,1)', loc, pars); 
    caxis([0 0.15])
    text(-0.5, 0.7,'mean coherence across frequencies (15-22 Hz)')
    
    %only for the strongest frequency
    subplot(2,2,2)
    showfield_general(l_COH(l_max_f, :)', loc, pars);
    text(-0.4, 0.7,'coherency at strongest frequency')
    
    %coherence spectrum at strongest channel
    subplot(2,2,3)
    plot(frqs(frq_inds), l_COH(:,l_max_chan))
    xlabel('Freq (Hz)')
    ylabel('Coherence')
    title('Coherence spectrum at strongest channel') 
    
    clear outname
    outname = sprintf('%scoh_analysis_Patient%s_left.png',DIROUT,patientID{id});
    print(outname,'-dpng');
    close all
    
    
    clear fileName D D_ft data conn COH absCOH loc l_COH r_COH
    
    
end