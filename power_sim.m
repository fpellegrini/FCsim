
clear all
% generation of sensor signal 
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
id = 1;

D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
load(sprintf('Filter_Patient%s.mat',patientID{id}));%Filter and whole CS
CS(:,:,47:end)=[];
id_meg_chan = 1:125;
id_meg_chan(D.badchannels)=[];
nmeg = numel(id_meg_chan);
inode = randi(size(A,2),1);

pv = fp_project_power(CS(id_meg_chan,id_meg_chan,:),A);
pow = squeeze(pv(inode,:));
clear L
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);

whitenoise = randn(size(pow));
whitenoise = whitenoise ./ norm(whitenoise, 'fro');
pow1 = 0.95*pow + 0.10*whitenoise; 
pow1 = pow1 ./ norm(pow1, 'fro');

L1 = squeeze(L(:,inode,1));

    
signal = L1*pow1;
    

%% project to source 

for ifreq= 1:46 
    signal_v(:,ifreq) = squeeze(A(:,:,ifreq))' * squeeze(signal(:,ifreq)); 
    signal_n(:,ifreq) = squeeze(A(:,:,ifreq))' * eye(size(signal(:,ifreq)));
end

a = signal_v; 

b = mean(a(:,2:10),2);
c = mean(a(:,11:17),2);
d = mean(a(:,2:17),2);

% b = squeeze(mean(mean(pow_noise,1),3));

outname = 'pow3.nii';
fp_data2nii(c,nan,[],outname)