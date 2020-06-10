function [mic, mim] = fp_get_mim_case2(A,CS,fqA,D)

%A is Beamformer filter
[nmeg, ni, nvox, nfreq] = size(A);
A2 = reshape(A,nmeg,ni*nvox,nfreq);

tic
CSroi = zeros(nvox*ni,nvox*ni,nfreq/4);
for ifq = 1: nfreq/4
    CSroi(:,:,ifq) = squeeze(A2(:,:,fqA(ifq)))' * CS(:,:,ifq)...
        * squeeze(A2(:,:,fqA(ifq)));
end
toc

%divide by power to obtain coherence
tic
Cohroi = zeros(size(CSroi));
for ifreq = 1: nfreq/4
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end
toc
clear CSroi


%% MIM only to aggegate 2 dimensions

npcs = repmat(2,[nvox,1]);
[mic_v,mim_v]= fp_mim(Cohroi,npcs);
clear Cohroi

%% sum up mim/ mic within rois

for iroi = 1:D.nroi 
    for jroi = 1:D.nroi
        mic(iroi,jroi,:) = squeeze(sum(sum(mic_v(D.ind_roi_cortex{iroi},...
            D.ind_roi_cortex{jroi},:),1),2));
        mim(iroi,jroi,:) = squeeze(sum(sum(mim_v(D.ind_roi_cortex{iroi},...
            D.ind_roi_cortex{jroi},:),1),2));
    end
end