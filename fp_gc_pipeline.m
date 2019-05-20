function fp_gc_pipeline 

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

if ~exist('shuffle','var')
    shuffle = 0;
end


for id = 1:numel(patientID)
    
    load(sprintf('Filter_Patient%s.mat',patientID{id}));%Filter and whole CS
    clear A 
    
    %include for shuffling and manual cs calc 
%     id_meg_trials = 1:n_trials;
%     
%     if shuffle == 1
%         rng('shuffle')
%         id_lfp_trials = randperm(n_trials);
%     else
%         id_lfp_trials = id_meg_trials;
%     end 
    
    
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    
    X = D(:,:,:);
    D_ft = ftraw(D);
    n_trials = length(D_ft.trial);
    
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frqs(frqs>90) = [];
    nfreq = numel(frqs);
    
    %construct filters 
    L1 = inverse.MEG.L;
    ns = numel(L1);
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    nmeg = numel(id_meg_chan);
    id_lfp_chan = 126:131;
    nlfp = numel(id_lfp_chan);
    for is=1:ns
        L(:,is,:)= L1{is};
    end    
    
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
   
    cCS = CS(1:(end-nlfp),end-nlfp+1:end,:);
    
    %project cross spectrum to voxel space
    for ifq = 1:nfreq
        for idir = 1:3
            CSv(idir,:,:,ifq) = squeeze(A(:,idir,:,ifq))' * cCS(:,:,ifq);
        end
    end
    
   
    
    
    
end 
    
    
    