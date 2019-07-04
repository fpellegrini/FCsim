clear all

id=5;

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

[~, voxID] = fp_find_commonvox;
nlags = 4;

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

%true CS
clear id_trials_1 id_trials_2 CS A_
id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(X,fres,'MT',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan], id_trials_1, id_trials_2);
    

%construct filters

load(sprintf('BF_Patient%s.mat',patientID{id}));
L1 = inverse.MEG.L;
ns_org = numel(L1);
for is=1:ns_org
    L(:,is,:)= L1{is};
end 


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

cCS = CS(1:(end-nlfp),end-nlfp+1:end,:); %nmeg x nlfp x nfreq
CSv = zeros(3*ns+nlfp,3*ns+nlfp,nfreq);
for ifq = 1:nfreq

    csv = zeros(ns*3+nlfp,ns*3+nlfp);
    csv(1:ns*3,end-nlfp+1:end) = squeeze(A_(:,:,ifq))' * cCS(:,:,ifq);
    csv(end-nlfp+1:end,1:ns*3)= csv(1:ns*3,end-nlfp+1:end)';
    csv(1:ns*3,1:ns*3) = squeeze(A_(:,:,ifq))' * CS(1:nmeg,1:nmeg,ifq) * squeeze(A_(:,:,ifq));
    csv(end-nlfp+1:end,end-nlfp+1:end) = CS(end-nlfp+1:end,end-nlfp+1:end,ifq);

    %replace power with real values
    clear n
    n = size(csv,1);
    csv(1:(n+1):end) = real(diag(csv));

    CSv(:,:,ifq) = csv; %.*10^4; %re-scale to avoid numerical errors
    clear csv
end

csvmeg = sum(real(CSv(1:end-nlfp,1:end-nlfp,:)),3); 
sfmeg = mean(diag(csvmeg));
csvlfp = sum(real(CSv(end-nlfp+1:end,end-nlfp+1:end,:)),3);
sflfp=mean(diag(csvlfp));

csvmeg = csvmeg./sfmeg;
csvlfp=csvlfp./sflfp;
%%

csvmeglfp =  CSv(1:end-nlfp,end-nlfp+1:end,:);
sfmeglfp = std(csvmeglfp(:));
csvmeglfp = csvmeglfp./sfmeglfp;

a=diag(csv);
%
%     %replace power with real values
%     clear n
%     n = size(csv,1);
%     csv(1:(n+1):end) = real(diag(csv));
%     CSv(:,:,ifq) = csv;
%     clear csv

figure
plot(squeeze(X(5,:,5)))
hold on 
plot(squeeze(X(126,:,5)))
legend('meg channel', 'lfp channel')

figure
imagesc(real(csv([3:5 end-2:end],[3:5 end-2:end])))
xlabel('3 meg chan + 3 lfp chan')
ylabel('3 meg chan + 3 lfp chan')
colorbar

%%

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

for id = 1: numel(patientID)
    
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    L1 = sources.L.MEG;
    ns_org = numel(L1);
    for is=1:ns_org
        L(:,is,:)= L1{is};
    end
    l = std(L(:));
    lr = range(L(:));
    clear L1 L
    
    
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    X = D(:,:,:);
    
    x(id)=std(X(:));
    clear X
    
end