function [sfmeg,sflfp] = fp_scaling_factor

DIROUT = '~/Dropbox/Franziska/Data_MEG_Project/';

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; % 

csvmeg=[];
csvlfp=[];
fres = 75;
csvmeg=zeros(125,1);
csvlfp=zeros(6,1);

for id = 1:numel(patientID)

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

    %true CS
    clear id_trials_1 id_trials_2 CS A_
    id_trials_1 = 1:n_trials;
    id_trials_2 = 1:n_trials;
    CS = fp_tsdata_to_cpsd(X,fres,'WELCH',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan], id_trials_1, id_trials_2);

    csvmeg(id_meg_chan) = csvmeg(id_meg_chan) + diag(sum(real(CS(1:end-nlfp,1:end-nlfp,:)),3)); 
    csvlfp = csvlfp + diag(sum(real(CS(end-nlfp+1:end,end-nlfp+1:end,:)),3));
  
end

sfmeg = sqrt(mean(csvmeg));
sflfp= sqrt(mean(csvlfp));

outname = sprintf('%sscaling_factor.mat',DIROUT);
save(outname,'sfmeg','sflfp','-v7.3')