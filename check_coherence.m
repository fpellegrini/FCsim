patientNumber = '04';

cd ~/Dropbox/MEG_Project/Data

fileName = '../../Data/redPLFP04_off.mat';

% spm_eeg_load(fileName)
D = wjn_meg_correct_mri(fileName);

D_ft = ftraw(D);

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

conn6 = data2spwctrgc(data, fres, nlags, cond, nboot, [], {'CS'});
cs6 = conn6.CS;
cs6 = permute(cs6,[2 3 1]);

for ifq = 1:46
    CSv6(ifq,:,:) = A(:,:,ifq)' * cs6(1:124,125:130,ifq);
end










   
    


