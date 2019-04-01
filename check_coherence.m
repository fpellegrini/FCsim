patientNumber = '04';

cd ~/Dropbox/MEG_Project/Data

fileName = '../../Data/redPLFP04_off.mat';

% spm_eeg_load(fileName)
D = wjn_meg_correct_mri(fileName);

D_ft = ftraw(D);

lfp_ind = 126;
meg_inds = 1:125;
meg_inds(D.badchannels)=[];
chan1= 100;
ifq = 5;

load(sprintf('Filter_Patient%s.mat','04'));

cs2 = CS(chan1,lfp_ind,ifq);
coh2 = fp_timesensor2sourcecoh(patientNumber, 0);
%%


N_trials = length(D_ft.trial);
[~, N_samples] = size(D_ft.trial{1});
N_chans = length(meg_inds);

fs = D.fsample;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90)=[];

nlags = 20; %does not influence the CS 
cond = 0; %does not influence the CS 
nboot = 1;

data = D(:,:,:);
data(132:end,:,:)=[];
data(D.badchannels,:,:) = [];

%%
data1 = data(1:5,:,:); 

conn6 = data2spwctrgc(data1, fres, nlags, cond, nboot, [], {'CS'});
cs6 = conn6.CS;

CSv6 = A(:,:,ifq)' * cs6(ifq,1:124,125)';
coh6 = cs2coh(CSv6);

coh2 = coh2(ifq,:,1)';








   
    


