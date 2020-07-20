function fp_fdr_mim(DIROUT)

alpha = 0.001;

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

%% load data

for id = 1:numel(patientID)
    load(sprintf('./roi_MIM90_sub%s.mat',patientID{id}))
    
    mim_true(id,:,:,:) = MIM_TRUE; 
    mim_shuf(id,:,:,:,:) = MIM_SHUF;
    mic_true(id,:,:,:) = MIC_TRUE; 
    mic_shuf(id,:,:,:,:) = MIC_SHUF;
    
    clear MIM_TRUE MIC_TRUE MIM_SHUF MIC_SHUF
end

[nsubs,nit,nroi,~,nfreq] = size(mim_shuf);


%% true

fprintf('Testing...\n')
[true_p_mic,~,~] = fp_get_signrank_results_megmeg(mic_true,mic_shuf,alpha);
[true_p_mim,~,~] = fp_get_signrank_results_megmeg(mim_true,mim_shuf,alpha);

[p_mic, mask_mic] = fdr( true_p_mic, 0.005);
[p_mim, mask_mim] = fdr( true_p_mim, 0.005);


%%


