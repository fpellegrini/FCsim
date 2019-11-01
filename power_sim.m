clear all

%parameters
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
id = 2;
nfreq = 46;
fres = 75;
inode = 2100; %randi(size(A,2),1);
isens = randi(125,1);
idir = randi(2,1);

%load MEG time series and CS
D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
X = D(:,:,:);
id_meg_chan = 1:125;
X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;
x = squeeze(X(isens,:,:)); %random time meg series
ntrials = size(x,2);
id_meg_trials = 1:ntrials;

%leadfield
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);
L1 = squeeze(L(:,inode,idir));
nmeg = size(L1,1);
id_meg_chan = 1:nmeg;

for itrial = 1:ntrials
    clear sig whitenoise 
    sig = L1 * x(:,itrial)';
    sig = sig ./ norm(sig, 'fro');
    %signal on sensor level 

    %add white noise 
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise; 
    signal(:,:,itrial) = sig ./ norm(sig, 'fro');
end

%filter
load(sprintf('Filter_Patient%s_e.mat',patientID{id}))
clear CS
ns = size(A,2);

%meg-meg CS
CS = fp_tsdata_to_cpsd(signal,fres,'MT',id_meg_chan, id_meg_chan, id_meg_trials, id_meg_trials);
CS(:,:,nfreq+1:end) = [];

%power
for ifreq = 1: nfreq   
    for is = 1:ns
        pow(is,ifreq) = real(squeeze(A(:,is,ifreq))' * CS(:,:,ifreq) * squeeze(A(:,is,ifreq)));
        pow_noise(is,ifreq) = real(squeeze(A(:,is,ifreq))' * eye(size(CS(:,:,ifreq))) * squeeze(A(:,is,ifreq)));
    end
end
%same as:
% pv = fp_project_power(CS,A);
    
%% project to source 

a = pow;% ./signal_n; 

b = mean(a(:,2:10),2);
c = mean(a(:,11:17),2);
d = mean(a(:,2:17),2);
e = mean(pow_noise,2);
f = mean(a,2);
g = mean(pow./pow_noise,2);
h = mean(pow - pow_noise*(pow_noise\pow),2);
subplot(1,2,1)
plot(f)
hold on 
plot(inode,0,'r+')
title('signal power')
xlabel('voxel id')
ylabel('pow')

subplot(1,2,2)
plot(e) 
title('noise power')
xlabel('voxel id')
ylabel('pow')

%%
outname = 'g_sub2_eloreta.nii';
fp_data2nii(g./10^-6,sources.pos,[],outname,id)