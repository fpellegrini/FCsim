function fp_sensorspace_coh(patientNumber)

cd ~/Dropbox/Data_MEG_Project/
DIROUT = '~/Dropbox/Data_MEG_Project/figures/sensorspace_coh/';
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

for id = [9:numel(patientID)] %patient 12: too few good channels for a plot
    
    fileName = sprintf('redPLFP%s_off', patientID{id});
    D = spm_eeg_load(fileName);
    D_ft = ftraw(D);
    N_trials = length(D_ft.trial);
    [~, N_samples] = size(D_ft.trial{1});
    meg_inds = MEG_INDS;
    meg_inds(D.badchannels) = [];
    N_meg = length(meg_inds);

    fs = D.fsample;
    fres = fs;
    frqs = sfreqs(fres, fs);
    
    data = [D_ft.trial{:}];
        
    data = reshape(data([meg_inds LFP_INDS], :), N_meg+N_lfp, N_samples, N_trials);

    conn = data2spwctrgc(data, fres, nlags, cond, nboot, [], {'CS'});        
    COH = cs2coh(conn.CS);
    absCOH = abs(COH(:, 1:N_meg, end-N_lfp+1:end));
    imCOH = abs(imag(COH(:, 1:numel(meg_inds), end-N_lfp+1:end)));
    
    loc = mk_sensors_plane(D_ft.grad.chanpos(meg_inds, [2 1 3]));
    frq_inds = find(frqs > 13 & frqs < 30);

    
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
    load(sprintf('perm_sensor_allchans_beta_Patient%s.mat',patientID{id}))

    %right lfp channels   
    r_COH = mean(absCOH(frq_inds,:,1:3),3); %15 x 125
    r_COH(:,p_r_abs>0.05)= nan;
    clear outname
    outname = sprintf('%scoh_abs_Patient%s_right.png',DIROUT,patientID{id});
    
    fp_plot_sensorspace_coh(r_COH,frqs, frq_inds, loc, [0 0.15] ,outname)
    
    %left lfp channels 
    l_COH = mean(absCOH(frq_inds,:,4:6),3); %15 x 125
    l_COH(:,p_l_abs>0.05)= nan;
    clear outname
    outname = sprintf('%scoh_abs_Patient%s_left.png',DIROUT,patientID{id});
    
    fp_plot_sensorspace_coh(l_COH,frqs,frq_inds, loc, [0 0.15],outname)
    
    %plot imaginary parts 
    %right
    rim_COH = mean(imCOH(frq_inds,:,1:3),3); %15 x 125
    rim_COH(:,p_r_im>0.05)= nan;
    clear outname
    outname = sprintf('%scoh_imag_Patient%s_right.png',DIROUT,patientID{id});
        
    fp_plot_sensorspace_coh(rim_COH,frqs, frq_inds, loc, [0 0.15],outname)
    
    %left
    lim_COH = mean(imCOH(frq_inds,:,4:6),3); %15 x 125
    lim_COH(:,p_l_im>0.05)= nan;
    clear outname
    outname = sprintf('%scoh_imag_Patient%s_left.png',DIROUT,patientID{id});
        
    fp_plot_sensorspace_coh(lim_COH,frqs, frq_inds, loc,[0 0.15],outname)
    
    clear fileName D D_ft data conn COH absCOH imCOH loc l_COH r_COH lim_COH rim_COH
    
    
end
end 


function fp_plot_sensorspace_coh(coh,frqs, frq_inds, loc,c_scale, outname)
    
    
    temp = max(max(coh));
    [max_f, max_chan] = find(coh==temp);
    clear temp
    
    figure
    figone(30,34)
    pars.scale = c_scale;
    %plot right channels
    %mean across freqs between 15 and 22
    subplot(2,2,1)
    showfield_general(mean(coh,1)', loc, pars); 
    caxis(c_scale)
    text(-0.5, 0.7,'Mean coherence across frequencies (15-22 Hz)')
    
    %only for the strongest frequency
    subplot(2,2,2)
    showfield_general(coh(max_f, :)', loc, pars);
    text(-0.4, 0.7,'Coherence at strongest frequency')
    
    %coherence spectrum at strongest channel
    subplot(2,2,3)
    plot(frqs(frq_inds), coh(:,max_chan))
    xlabel('Freq (Hz)')
    ylabel('Coherence')
    title('Coherence spectrum at strongest channel') 
    
    print(outname,'-dpng');
    close all
end 
    