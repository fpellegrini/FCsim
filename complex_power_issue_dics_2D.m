
%% short version

clear all
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[commonvox_pos, voxID] = fp_find_commonvox;

for id = 1:numel(patientID)
    clearvars -except id patientID pow voxID pow_noise commonvox_pos
    
    load(sprintf('Filter_Patient%s_dics_2D_wo_reg.mat',patientID{id}))  
    
    mni_pos = fp_getMNIpos(patientID{id});
    [~, noEq] = fp_symmetric_vol(mni_pos);
    A(:,:,noEq,:) = [];
    filter = A(:,:,voxID{id},:);
    
    [ndim, nmeg,ns,nfreq] = size(filter);
    CS = CS(1:nmeg,1:nmeg,:);
    %power
    for ifreq = 1:nfreq
        for idim = 1:ndim
            for is = 1:ns
                pow(id,is,idim,ifreq) = real(squeeze(A(idim,:,is,ifreq)) * CS(:,:,ifreq) * squeeze(A(idim,:,is,ifreq))');
            end
        end
    end
end

e = squeeze(sum(sum(sum(pow,1),3),4));
beta = squeeze(sum(sum(sum(pow(:,:,:,6:15),1),3),4));
%%
outname = 'real_pow_beta_dics2D.nii';
fp_data2nii(beta./10^-8,commonvox_pos,[],outname)

outname = 'real_pow_all_dics2D.nii';
fp_data2nii(e./10^-8,commonvox_pos,[],outname)

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
% CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
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