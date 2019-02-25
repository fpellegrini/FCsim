cd ~/Dropbox/MEG_Project/Data
load('BF_Patient08.mat')
X = data.D(:,:,:);
fres=300;
id_meg_chan = 1:125;
id_lfp_chan = 126:131;
id_meg_trials = 1:size(X,3);
id_lfp_trials = id_meg_trials; %randperm(size(X,3)); %

S_org = fp_tsdata_to_cpsd(X,[],fres,id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials);


id_lfp_trials = randperm(size(X,3)); %
S_shuf = fp_tsdata_to_cpsd(X,S_org,fres,id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials);

figure
imagesc(squeeze(abs(S_org(id_meg_chan,id_lfp_chan,10))))
title('original')

figure
imagesc(squeeze(abs(S_shuf(id_meg_chan,id_lfp_chan,10))))
title('shuffled')

figure
imagesc(squeeze(a(1:125,126:131,11)))