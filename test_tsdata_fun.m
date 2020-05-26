
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
fres = 50; 
for id = 1:numel(patientID)
    
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    X = D(id_meg_chan,:,:);
    ndat = size(X,2);
    
    tic
    clear CS tmp
    CS = tsdata_to_cpsd(X, fres, 'MT',ndat);
    for ifreq = 1:fres+1 
        tmp(:,ifreq) = diag(CS(:,:,ifreq));
    end
    power_mt_nw3(id,:)=squeeze(mean(tmp,1));
    toc
    
    tic
    clear CS tmp
    CS = tsdata_to_cpsd(X, fres, 'WELCH',ndat);
    for ifreq = 1:fres+1 
        tmp(:,ifreq) = diag(CS(:,:,ifreq));
    end
    power_welch(id,:) = squeeze(mean(tmp,1)); 
    toc
    
    tic
    clear CS tmp
    CS = tsdata_to_cpsd(X,fres,'MT',ndat,[],1);
    for ifreq = 1:fres+1 
        tmp(:,ifreq) = diag(CS(:,:,ifreq));
    end
    power_mt_nw1(id,:) = squeeze(mean(tmp,1)); 
    toc
    
    tic
    clear tmp
    for itrial = 1:size(X,3)
        [tmp(:,:,itrial), f] = pwelch(X(:,:,itrial)',ndat, ndat/2, 2*fres,D.fsample);
    end
    power_pwelch(id,:) = squeeze(mean(mean(tmp,3),2));
    toc
end
    