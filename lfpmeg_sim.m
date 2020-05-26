clear all

%% signal generation
%parameters

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2;
% nfreq = 46;
delay = 30; %in samples, is equal to 100 ms
inode = 2100; %randi(size(A,2),1);
idir = randi(2,1);
lfp_chan = 127;
nmeg = 125; 

%load MEG time series and CS
D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
X = D(:,:,:);
id_meg_chan = 1:125;
X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;
fs = D.fsample; 
 
%correlated nodes 
x1 = squeeze(X(lfp_chan,:,:)); %random time meg series
n_trials = size(x1,2);
x2(1,:,:) = x1(delay+1:end,:); %time series with time delay

%filter signal in alpha band
[b, a] = butter(2, [5 12]/(fs/2),'pass');
x = squeeze(reshape(filtfilt(b, a, reshape(x2, 1, [])')',1,[],n_trials));
x1=x;

%uncorrelated nodes
% for ii =1:n_rand_nodes
%     ind = randperm(70);
%     xx(ii+2,:,:) = squeeze(X(isens,delay:end,ind));
% end

%leadfield
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);
L1 = squeeze(L(:,[inode],idir));
id_meg_chan = 1:nmeg;

%project to sensors and add white noise 
for itrial = 1:n_trials
    clear sig whitenoise
    sig = L1 * x1(:,itrial)';
    sig = sig ./ norm(sig, 'fro');
    %signal on sensor level
    
    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal(:,:,itrial) = sig ./ norm(sig, 'fro');
end

signal = cat(1,signal,X(nmeg+1:end,1:end-delay,:));

clear x1 x2 xx whitenoise sig X 
%% lfp meg pipeline 

fres = 75;
id_lfp_chan = 126:131;
id_meg_trials = 1:n_trials;
id_lfp_trials = id_meg_trials;
nlfp = numel(id_lfp_chan);
ns = size(L,2);
nfreq = 46;
[commonvox_pos, voxID] = fp_find_commonvox;

%CS
CS = fp_tsdata_to_cpsd(signal,fres,'WELCH',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan],...
    id_meg_trials, id_lfp_trials);
CS(:,:,nfreq+1:end)=[];

%filter
A=zeros(nmeg,ns,nfreq);

for ifrq = 1:nfreq
    A(:,:,ifrq) = fp_filter(squeeze(CS(1:end-nlfp,1:end-nlfp,ifrq)), L);
end

cCS = CS(1:(end-nlfp),end-nlfp+1:end,:);

%project cross spectrum to voxel space
for ifq = 1:nfreq
    CSv(ifq,:,:) = A(:,:,ifq)' * cCS(:,:,ifq);
end

%get voxel power
pv = fp_project_power(CS(1:nmeg,1:nmeg,:),A);

%coherence
coh = CSv;
for ifreq = 1:nfreq
    clear plfp
    plfp = diag(squeeze(CS(nmeg+1:end,nmeg+1:end,ifreq))); %power of lfp channels
    coh(ifreq, :, :) = squeeze(CSv(ifreq, :, :)) ...
        ./ sqrt(pv(:,ifreq)*plfp');
end

mni_pos = fp_getMNIpos(patientID{id});

%get flip id and symmetric head
[sym_pos, noEq] = fp_symmetric_vol(mni_pos);
match_pos = sym_pos(voxID{id},:);
[~,flip_id] = fp_flip_vol(match_pos);

coh(:,noEq,:) = [];
match_coh = coh(:,voxID{id},:);
flip_coh = match_coh;
flip_coh(:,:,4:6) = coh(:,flip_id,4:6);

abs_coh = abs(imag(flip_coh));

%%
a = squeeze(mean(abs_coh,3));

a1 = squeeze(mean(a(3:6,:),1));
a2 = squeeze(mean(a(10:end,:),1));

outname = sprintf('1.nii');
% outname = sprintf('lfpmeg_sim_dics_%d.nii',inode);
fp_data2nii(a1,commonvox_pos,[],outname,id)

outname = sprintf('2.nii');
fp_data2nii(a2,commonvox_pos,[],outname,id)




