function [mim, mic, mean_coh,to_save] = fp_correct_mim(A,signal_sensor_save, fqA, nfqA, D, ihemi, mic, mim, mean_coh, to_save)

nit = 50; 
n_trials = size(signal_sensor_save,3); 


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
    
    clear mim1 mic1 mean_coh1
    [mic_max_shuf(:,:,:,iit), mim_max_shuf(:,:,:,iit), ~, mean_coh_max_shuf(:,:,:,iit)] = fp_get_mim(A,CS,fqA,nfqA, D,ihemi,'max');
    
    clear mim1 mic1 mean_coh1
    [mic90_shuf(:,:,:,iit), mim90_shuf(:,:,:,iit), ~, mean_coh90_shuf(:,:,:,iit)] = fp_get_mim(A,CS,fqA,nfqA, D,ihemi,'percent');
    
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
mean_coh.max_corrected = (mean_coh.max - m_max)./s_max;

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
mean_coh.percent_corrected = (mean_coh.percent - m90)./s90; 

to_save.max.mim_shuf = mim_max_shuf; 
to_save.max.mic_shuf = mic_max_shuf; 
to_save.max.mean_coh_shuf = mean_coh_shuf; 

to_save.percent.mim_shuf = mim_percent_shuf; 
to_save.percent.mic_shuf = mic_percent_shuf; 
to_save.percent.mean_coh_shuf = mean_coh_shuf; 




