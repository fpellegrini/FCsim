function [mic, mim, to_save, mean_icoh,mean_acoh] = fp_correct_mim(A,signal_sensor, fqA, nfqA, D, ihemi, mic, mim, mean_icoh, mean_acoh, to_save,fres)

nit = 50; 
n_trials = size(signal_sensor,3); 
zs=1;


for iit = 1:nit 
    
    fprintf('Iteration %d... \n',iit)
    id_trials_1 = 1:n_trials;
    id_trials_2 = randperm(n_trials);
    id_meg_chan = 1:size(signal_sensor,1);
    nmeg = numel(id_meg_chan);
    
    %cross spectrum
    fprintf('Calculating cross spectrum... \n')
    CS = fp_tsdata_to_cpsd(signal_sensor,fres,'WELCH',...
        [id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
    CS(:,:,1)=[];
    nfreq = size(CS,3);
    
    [mic_max_shuf(:,:,:,iit), mim_max_shuf(:,:,:,iit), ~, ...
        mean_coh_max_shuf(:,:,:,iit),mean_abscoh_max_shuf(:,:,:,iit)] = ...
        fp_get_mim(A,CS,fqA,nfqA, D,ihemi,'max',zs);
    
    [mic90_shuf(:,:,:,iit), mim90_shuf(:,:,:,iit), ~, ...
        mean_coh90_shuf(:,:,:,iit), mean_abscoh90_shuf(:,:,:,iit) ] = ...
        fp_get_mim(A,CS,fqA,nfqA, D,ihemi,'percent',zs);
    
end

%%

clear m_max s_max
m_max = mean(mic_max_shuf,4); 
s_max = std(mic_max_shuf,0,4); 
mic.max_corrected = (mic.max - m_max)./s_max;

clear m_max s_max
m_max = mean(mim_max_shuf,4); 
s_max = std(mim_max_shuf,0,4); 
mim.max_corrected = (mim.max - m_max)./s_max;

clear m_max s_max
m_max = mean(mean_coh_max_shuf,4); 
s_max = std(mean_coh_max_shuf,0,4); 
mean_icoh.max_corrected = (mean_icoh.max - m_max)./s_max;

clear m_max s_max
m_max = mean(mean_abscoh_max_shuf,4); 
s_max = std(mean_abscoh_max_shuf,0,4); 
mean_acoh.max_corrected = (mean_acoh.max - m_max)./s_max;

clear m90 s90
m90 = mean(mic90_shuf,4); 
s90 = std(mic90_shuf,0,4); 
mic.percent_corrected = (mic.percent - m90)./s90; 

clear m90 s90
m90 = mean(mim90_shuf,4); 
s90 = std(mim90_shuf,0,4); 
mim.percent_corrected = (mim.percent - m90)./s90; 

clear m90 s90
m90 = mean(mean_coh90_shuf,4); 
s90 = std(mean_coh90_shuf,0,4); 
mean_icoh.percent_corrected = (mean_icoh.percent - m90)./s90; 

clear m90 s90
m90 = mean(mean_abscoh90_shuf,4); 
s90 = std(mean_abscoh90_shuf,0,4); 
mean_acoh.percent_corrected = (mean_acoh.percent - m90)./s90; 

to_save.max_corrected.mim_shuf = mim_max_shuf; 
to_save.max_corrected.mic_shuf = mic_max_shuf; 
to_save.max_corrected.mean_icoh_shuf = mean_coh_max_shuf; 
to_save.max_corrected.mean_acoh_shuf = mean_abscoh_max_shuf; 

to_save.percent_corrected.mim_shuf = mim90_shuf; 
to_save.percent_corrected.mic_shuf = mic90_shuf; 
to_save.percent_corrected.mean_icoh_shuf = mean_coh90_shuf; 
to_save.percent_corrected.mean_acoh_shuf = mean_abscoh90_shuf; 

%%
nvoxroi_all = to_save.nvoxroi'* to_save.nvoxroi;
nvoxroi_all = nvoxroi_all(:);

c1 = sum(mim.max_corrected,3);
c2 = sum(mic.max_corrected,3);
c3 = sum(mean_icoh.max_corrected,3);
c4 = sum(mean_acoh.max_corrected,3);
to_save.max_corrected.corr_voxmim = corr(nvoxroi_all,c1(:));
to_save.max_corrected.corr_voxmic = corr(nvoxroi_all ,c2(:));
to_save.max_corrected.corr_voxnpcs = corr(to_save.nvoxroi', to_save.max.npcs');
to_save.max_corrected.corr_voxmeancoh = corr(nvoxroi_all,c3(:));
to_save.max_corrected.corr_voxmeanabscoh = corr(nvoxroi_all,c4(:));

c1 = sum(mim.percent_corrected,3);
c2 = sum(mic.percent_corrected,3);
c3 = sum(mean_icoh.percent_corrected,3);
c4 = sum(mean_acoh.percent_corrected,3);
to_save.percent_corrected.corr_voxmim = corr(nvoxroi_all,c1(:));
to_save.percent_corrected.corr_voxmic = corr(nvoxroi_all,c2(:));
to_save.percent_corrected.corr_voxnpcs = corr(to_save.nvoxroi', to_save.percent.npcs');
to_save.percent_corrected.corr_voxmeancoh = corr(nvoxroi_all,c3(:));
to_save.percent_corrected.corr_voxmeanabscoh = corr(nvoxroi_all,c4(:));





