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

%meg-meg CS
CS = fp_tsdata_to_cpsd(signal,fres,'MT',id_meg_chan, id_meg_chan, id_meg_trials, id_meg_trials);
CS(:,:,nfreq+1:end) = [];

ns = size(L,2);

%filter
A=zeros(2,nmeg,ns,nfreq);

for ifrq = 1:nfreq
    cCS = CS(:,:,ifrq);
    lambda = mean(diag(real(cCS)))/100;
    
    CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
    
    for is=1:ns %iterate across nodes
        Lloc=squeeze(L(:,is,:));
        A(:,:,is,ifrq) = pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv; %create filter
    end
end


%power
for idim = 1:2
    for ifreq = 1: nfreq
        for is = 1:ns
            pow(is,idim,ifreq) = real(squeeze(A(idim,:,is,ifreq)) * CS(:,:,ifreq) * squeeze(A(idim,:,is,ifreq))');
            pow_noise(is,idim,ifreq) = real(squeeze(A(idim,:,is,ifreq)) * eye(size(CS(:,:,ifreq))) *  squeeze(A(idim,:,is,ifreq))');
        end
    end
end

pow = squeeze(sum(pow,2));
pow_noise = squeeze(sum(pow_noise,2));
    
%% project to source 

g = pow ./ pow_noise;
g = squeeze(mean(g,2));

plot(g)
hold on 
plot(inode,0,'r+')
title('power')
xlabel('voxel id')
ylabel('pow')

%%
outname = 'g_sub2_dics2D.nii';
fp_data2nii(g./10^-5,sources.pos,[],outname,id)