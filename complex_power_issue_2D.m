
%% short version

clear all
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[commonvox_pos, voxID] = fp_find_commonvox;

for id = 1:numel(patientID)
    clearvars -except id patientID pow voxID pow_noise commonvox_pos
    
    load(sprintf('Filter_Patient%s_e2D.mat',patientID{id}))  
    
    mni_pos = fp_getMNIpos(patientID{id});
    [~, noEq] = fp_symmetric_vol(mni_pos);
    A(:,noEq,:) = [];
    filter = A(:,voxID{id},:);
    
    [nmeg, ns, ndim] = size(filter);
    CS = CS(1:nmeg,1:nmeg,:);
    %power
    for ifreq = 1:size(CS,3)
        for idim = 1:2
            for is = 1:ns
                pow(id,is,idim,ifreq) = real(squeeze(A(:,is,idim))' * CS(:,:,ifreq) * squeeze(A(:,is,idim)));
            end
        end
    end
end

e = squeeze(sum(sum(sum(pow,1),3),4));
beta = squeeze(sum(sum(sum(pow(:,:,:,6:15),1),3),4));
%%
outname = 'real_pow_beta_e2D.nii';
fp_data2nii(beta./10^-8,commonvox_pos,[],outname)

outname = 'real_pow_all_e2D.nii';
fp_data2nii(e./10^-8,commonvox_pos,[],outname)

