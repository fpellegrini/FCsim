function fp_check_lfp_chans

cd ~/Dropbox/Data_MEG_Project/
DIROUT = '~/Dropbox/Data_MEG_Project/figures/lfp_chans/';
if ~exist(DIROUT); mkdir(DIROUT); end 

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors
fs = 300;

for id = 1:numel(patientID)
    
    fileName = sprintf('redPLFP%s_off', patientID{id});
    D = spm_eeg_load(fileName);
    D_ft = ftraw(D);
    N_trials = length(D_ft.trial);
    
    clear Pxx
    for itrial = 1:N_trials
        [Pxx(:, :, itrial), F] = pwelch(zscore(D_ft.trial{itrial}(126:131, :)'), [], [], 4*fs, fs);
    end
    
    figure
    semilogy(F, mean(Pxx, 3)); grid on %mean across trials
    legend(D.chanlabels{126},D.chanlabels{127},D.chanlabels{128},...
        D.chanlabels{129},D.chanlabels{130},D.chanlabels{131})
    ylim([10^(-6) 10^0])
    xlabel('Freq. [Hz]')
    ylabel('PSD [dB]')
    
    outname = sprintf('%s%s.png',DIROUT,patientID{id});
    print(outname,'-dpng');
    close all
    
    
end