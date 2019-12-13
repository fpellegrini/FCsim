
DIRDATA = './';
DIRFIG = '~/Dropbox/Franziska/Data_MEG_Project/alpha_sensorlevel/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};


for id = [1:6 8:numel(patientID)] 
    clearvars -except DIRDATA patientID id DIRFIG
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    fs = D.fsample;
    
    X = D(:,:,:);
    D_ft = ftraw(D);
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;

    clear Pxx
    for itrial = 1: size(X,3)
        [Pxx(:,:,itrial), f] = pwelch(X(id_meg_chan,:,itrial)', [], [], 4*fs, fs);
    end
    
    a =mean(mean(Pxx,3),2);
    b = mean(mean(Pxx(33:49,:,:),1),3);
    c = mean(mean(Pxx(:,:,:),1),3);
    
    figure
    semilogy(f, a); grid on %mean across occipital channels
    xlabel('Freq. [Hz]')
    ylabel('PSD [dB]')
    
    outname = sprintf('%sspectrum_sub_%s',DIRFIG,patientID{id});
    print(outname,'-dpng')
    close all
    
    chanpos = D_ft.grad.chanpos;
    chanpos(D.badchannels,:)=[];
    
    loc = mk_sensors_plane(chanpos(:,[2 1 3]));

    figure
    showfield_general(log(b),loc);
    colormap('jet')
    
    outname = sprintf('%salphatopo_sub_%s',DIRFIG,patientID{id});
    print(outname,'-dpng')
    close all
    
    
    figure
    showfield_general(log(c),loc);
    colormap('jet')
    
    outname = sprintf('%salltopo_sub_%s',DIRFIG,patientID{id});
    print(outname,'-dpng')
    close all
    
    %%
end
