function [mic_2, mim_2, baseline] = fp_get_mim_case2(A,CS,fqA,D)
% case two and benchmark mim

%A is Beamformer filter
[nmeg, ni, nvox, nfreq] = size(A);
A2 = reshape(A,nmeg,ni*nvox,nfreq);

%CS on voxel level 
tic
CSroi = zeros(nvox*ni,nvox*ni,nfreq);
for ifq = 1: nfreq
    CSroi(:,:,ifq) = squeeze(A2(:,:,fqA(ifq)))' * CS(:,:,ifq)...
        * squeeze(A2(:,:,fqA(ifq)));
end
toc

%divide by power to obtain coherence
tic
Cohroi = zeros(size(CSroi));
for ifreq = 1: nfreq
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end
toc
clear CSroi

%% baseline
%in sub_ind_roi, there are the 'true' voxels of the signal generation 

baseline = Cohroi(D.sub_ind_cortex,D.sub_ind_cortex,:);
iReg = length(D.sub_ind_roi{1});
if iReg ~=1 %mim across two voxels of one region 
    [baseline_mic,baseline_mim]= fp_mim(baseline,repmat(iReg,D.nroi,1));
    clear baseline 
    baseline.mic = baseline_mic; 
    baseline.mim = baseline_mim; 
    clear baseline_mic baseline_mim
end



%% case two: MIM only to aggegate 2 dimensions

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