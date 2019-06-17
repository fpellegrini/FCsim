clear all

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id=10;
[~, voxID] = fp_find_commonvox;
nlags = 4;
load(sprintf('Filter_Patient%s.mat',patientID{id}));%load true CS for constructing filters
clear A

%load data
D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
X = D(:,:,:);
D_ft = ftraw(D);
n_trials = length(D_ft.trial);

%channel IDs
id_meg_chan = 1:125;
id_meg_chan(D.badchannels)=[];
nmeg = numel(id_meg_chan);
id_lfp_chan = 126:131;
nlfp = numel(id_lfp_chan);

%frequency parameters
fs = D.fsample;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>10) = [];
nfreq = numel(frqs);
freqs = linspace(0, 1, fres+1);
maxfreq = fres+1;
freqs = freqs(1:maxfreq);
z = exp(-i*pi*freqs);

%construct filters

load(sprintf('BF_Patient%s.mat',patientID{id}));
L1 = inverse.MEG.L;
ns_org = numel(L1);
for is=1:ns_org
    L(:,is,:)= L1{is};
end 
L = L.* (10^(-log10(range(L(:)))+5.5));

%delete voxels that are not common in all subs
mni_pos = fp_getMNIpos(patientID{id});
[~, noEq] = fp_symmetric_vol(mni_pos);
L(:,noEq,:) = [];
L = L(:,voxID{id},:);
ns = numel(voxID{id});

A = nan(nmeg,3,ns,nfreq);
for ifrq = 1:nfreq
    clear currentCS lambda CSinv
    currentCS = squeeze(CS(1:end-nlfp,1:end-nlfp,ifrq));
    lambda = mean(diag(real(currentCS)))/100;
    CSinv=pinv(real(currentCS)+lambda * eye(size(currentCS)));
    
    for is=1:ns %iterate across nodes
        clear Lloc
        Lloc=squeeze(L(:,is,:));
        A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
    end
end

A_ = reshape(A, [nmeg, 3*ns, nfreq]);
clear A
iit=1;

cCS = CS(1:(end-nlfp),end-nlfp+1:end,:); %nmeg x nlfp x nfreq
CSv = zeros(3*ns+nlfp,3*ns+nlfp,nfreq);
ifq = 3;

csv = zeros(ns*3+nlfp,ns*3+nlfp);
csv(1:ns*3,end-nlfp+1:end) = squeeze(A_(:,:,ifq))' * cCS(:,:,ifq);
csv(end-nlfp+1:end,1:ns*3)= csv(1:ns*3,end-nlfp+1:end)';
csv(1:ns*3,1:ns*3) = squeeze(A_(:,:,ifq))' * CS(1:nmeg,1:nmeg,ifq) * squeeze(A_(:,:,ifq));
csv(end-nlfp+1:end,end-nlfp+1:end) = CS(end-nlfp+1:end,end-nlfp+1:end,ifq);

csv=csv.*10^8;

a=diag(csv);
%
%     %replace power with real values
%     clear n
%     n = size(csv,1);
%     csv(1:(n+1):end) = real(diag(csv));
%     CSv(:,:,ifq) = csv;
%     clear csv



imagesc(real(csv(1:3,[1:3 end-3:end])))
colorbar