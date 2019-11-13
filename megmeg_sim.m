clear all

%% signal generation 
%parameters
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
id = 2;
% nfreq = 46;
fres = 75;
delay = 30; %in samples, is equal to 100 ms 
inode1 = 2100; %randi(size(A,2),1);
inode2 = 1000;
isens = randi(125,1);
idir = randi(2,1);

%load MEG time series and CS
D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
X = D(:,:,:);
id_meg_chan = 1:125;
X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;

x1 = squeeze(X(isens,:,:)); %random time meg series
x2 = x1(delay:end,:); %second time series with time delay
x1 = x1(1:end-delay+1,:);

xx(1,:,:)= x1;
xx(2,:,:) = x2; 
n_trials = size(x1,2);
id_meg_trials = 1:n_trials;

%leadfield
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);
L1 = squeeze(L(:,[inode1 inode2],idir));
nmeg = size(L1,1);
id_meg_chan = 1:nmeg;

for itrial = 1:n_trials
    clear sig whitenoise 
    sig = L1 * xx(:,:,itrial);
    sig = sig ./ norm(sig, 'fro');
    %signal on sensor level 

    %add white noise 
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise; 
    signal(:,:,itrial) = sig ./ norm(sig, 'fro');
end

clear x1 x2 xx whitenoise sig X L L1

%% megmeg pipeline start 

id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(signal,fres,'MT',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);

%leadfield
clear L
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);
ns_org = size(L,2);

clear mni_pos label code roi_id u_roi_id csroi
mni_pos = fp_getMNIpos(patientID{id});
for ii = 1: ns_org
    [~,~,roi_id(ii)]=fp_get_mni_anatomy_new(mni_pos(ii,:));
end
u_roi_id = sort(unique(roi_id));
nroi = numel(u_roi_id)-1; %because white voxels are not counted





