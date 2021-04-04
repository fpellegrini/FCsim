
%% short version

clear all
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[commonvox_pos, voxID] = fp_find_commonvox;
pow = nan(numel(patientID),numel(voxID{1}),45);
pow_noise = nan(numel(patientID),numel(voxID{1}),45);
fres = 75;


for id = 1:numel(patientID)
    clearvars -except id patientID pow voxID pow_noise commonvox_pos pow_eloreta fres
    
    load(sprintf('Filter_Patient%s_e.mat',patientID{id})) 
    clear CS 
    
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    X = D(:,:,:);
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;
    id_trials = 1:size(X,3);
    nfreq = size(A,3);
    
    CS = fp_tsdata_to_cpsd(X,fres,'WELCH',id_meg_chan, id_meg_chan, id_trials, id_trials);
    CS(:,:,nfreq+1:end) = [];
    
    mni_pos = fp_getMNIpos(patientID{id});
    [~, noEq] = fp_symmetric_vol(mni_pos);
    A(:,noEq,:) = [];
    filter = A(:,voxID{id},:);
    filter(:,:,1)=[];
    
    [nmeg, ns, nfreq] = size(filter);  

    for ifreq = 1: nfreq    

        cfilter = filter(:,:,ifreq)';

        for is = 1:ns    
            pow(id,is,ifreq) = real(cfilter(is,:) * CS(:,:,ifreq) * cfilter(is,:)'); 
            pow_noise(id,is,ifreq) = real(cfilter(is,:) * (eye(size(CS(:,:,ifreq))).*1) * cfilter(is,:)'); 
        end    
    end 
end

%%
% pow = pow./pow_noise;

% %%
% ae = squeeze(sum(pow_eloreta,1));
% ad = squeeze(sum(pow,1));
% sae = sum(ae(:,4:6),2);
% sad = sum(ad(:,4:6),2);
% 
% %%
% 
% scatter3(commonvox_pos(:,1),commonvox_pos(:,2),commonvox_pos(:,3),10,A__)
% figure
% scatter3(commonvox_pos(:,1),commonvox_pos(:,2),commonvox_pos(:,3),10,As__)

h = squeeze(mean(pow(:,:,4:6),3));
% h1 = squeeze(mean(h,1));

for id = 1:numel(patientID)
    outname = sprintf('real_pow_alpha_eloreta_patient%s.nii',patientID{id});
    fp_data2nii(squeeze(h(id,:))./std(h(id,:)),commonvox_pos,[],outname,[])
end

% outname = 'pow2_e.nii';
% fp_data2nii(c,nan,[],outname)
% 
% outname = 'pow3_e.nii';
% fp_data2nii(d,nan,[],outname)


%% version from scratch, with 3D filters (outdated!) 

% D = spm_eeg_load('redPLFP10_off');
% 
% X = D(:,:,:);
% X = X./10^(log10(range(X(:)))-2);
% 
% id_meg_chan = 1:125;
% id_meg_chan(D.badchannels)=[];
% nmeg = numel(id_meg_chan);
% 
% fs = D.fsample;
% fres = 75;
% frqs = sfreqs(fres, fs);
% frqs(frqs>90) = [];
% nfreq = numel(frqs);
% 
% id_trials_1 = 1: size(X,3);
% id_trials_2 = id_trials_1;
% CS = fp_tsdata_to_cpsd(X,fres,'WELCH',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
% 
% %%%%%up to this point, the power still consists of doubles
% 
% %construct filters
% 
% %leadfield
% load('BF_Patient10.mat');
% L1 = inverse.MEG.L;
% ns = numel(L1);
% for is=1:ns
%     L(:,is,:)= L1{is};
% end
% L = L.* (10^(-log10(range(L(:)))));
% 
% A = nan(nmeg,ns,nfreq);
% for ifrq = 1:nfreq
%     clear currentCS lambda CSinv
%     currentCS = squeeze(CS(:,:,ifrq)); %nmeg x nmeg x nfq
%     lambda = mean(diag(real(currentCS)))/100;
%     CSinv=pinv(real(currentCS)+lambda * eye(size(currentCS)));
%     
%     for is=1:ns %iterate across nodes
%         clear Lloc
%         Lloc=squeeze(L(:,is,:));
%         filter = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv); %create filter
%         
%         csd = filter*real(squeeze(CS(:,:,ifrq)))*filter';
%         [u,~,~] = svd(csd);
%         LF = Lloc*u(:,1);
% 
%         %recompute filter in best orientation 
%         A(:,is,ifrq)=pinv((LF'*CSinv*LF))*LF'*CSinv;
%     end
% end
% 
% 
% %project cross spectrum to voxel space and get power and coherence
% for ifq =1:nfreq
%         
%         CSv = squeeze(A(:,:,ifq))' * ...
%             CS(:,:,ifq) * squeeze(A(:,:,ifq)); %3dim x nvox x nvox x nfreq
%         
%        pv(:,ifq) = real(diag(squeeze(CSv)));
%         
%         
% end
% 
% betapow = mean(pv(:,7:16),2);
% % fp_data2nii(betapow,1, [], 'betapow04.nii')