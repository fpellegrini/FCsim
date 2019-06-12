
%% short version
load('Filter_Patient18.mat') %load filter and true CS for subject 18
[nmeg, ns, nfreq] = size(A);
CS = CS(1:nmeg,1:nmeg,:);

pow = nan(ns,nfreq);

for ifreq = 2:nfreq

    cfilter = A(:,:,ifreq)';
    
    for is = 1:ns
        pow(is,ifreq) = cfilter(is,:) * CS(:,:,ifreq) * cfilter(is,:)';
    end
end

%% version from scratch, with 3D filters

D = spm_eeg_load('redPLFP10_off');

X = D(:,:,:);
X = X./10^(log10(range(X(:)))-2);

id_meg_chan = 1:125;
id_meg_chan(D.badchannels)=[];
nmeg = numel(id_meg_chan);

fs = D.fsample;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90) = [];
nfreq = numel(frqs);

id_trials_1 = 1: size(X,3);
id_trials_2 = id_trials_1;
CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);

%%%%%up to this point, the power still consists of doubles

%construct filters

%leadfield
load('BF_Patient10.mat');
L1 = inverse.MEG.L;
ns = numel(L1);
for is=1:ns
    L(:,is,:)= L1{is};
end
L = L.* (10^(-log10(range(L(:)))));

A = nan(nmeg,ns,nfreq);
for ifrq = 1:nfreq
    clear currentCS lambda CSinv
    currentCS = squeeze(CS(:,:,ifrq)); %nmeg x nmeg x nfq
    lambda = mean(diag(real(currentCS)))/100;
    CSinv=pinv(real(currentCS)+lambda * eye(size(currentCS)));
    
    for is=1:ns %iterate across nodes
        clear Lloc
        Lloc=squeeze(L(:,is,:));
        filter = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv); %create filter
        
        csd = filter*real(squeeze(CS(:,:,ifrq)))*filter';
        [u,~,~] = svd(csd);
        LF = Lloc*u(:,1);

        %recompute filter in best orientation 
        A(:,is,ifrq)=pinv((LF'*CSinv*LF))*LF'*CSinv;
    end
end


%project cross spectrum to voxel space and get power and coherence
for ifq =1:nfreq
        
        CSv = squeeze(A(:,:,ifq))' * ...
            CS(:,:,ifq) * squeeze(A(:,:,ifq)); %3dim x nvox x nvox x nfreq
        
       pv(:,ifq) = real(diag(squeeze(CSv)));
        
        
end

betapow = mean(pv(:,7:16),2);
% fp_data2nii(betapow,1, [], 'betapow04.nii')