function fp_cluster_mim(DIROUT,testmethod,alpha,fwf,j)

%Group statistics on megmeg data.
%Clustering group statistics across space and frequencies.
%
%Input:
%DIROUT
%
%testmethod = either 's' for signrank testing or 't' for simple
%thresholding
%
%fwf: testing method (less conservative): test first cluster against
%first, second with second etc. By default, clusters are always compared
%against the largest shuffled cluster.
%
%j=1 : only Julian's subjects


fp_addpath

if ~exist(DIROUT); mkdir(DIROUT); end

alpha_s = num2str(alpha);
alpha_s(1:2)=[];

if fwf==0
    fwf_s = [];
elseif fwf ==1
    fwf_s = 'fwf';
else
    error('Wrong fwf input')
end

if j == 0
    j_s = 'allsubs';
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
elseif j == 1
    j_s = 'j';
    patientID = {'04'; '07'; '09'; '10';'11';'20';'22';'25'};
else
    error('Wrong input j')
end

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
if strcmp(testmethod,'s')
    tic
    [true_p_mic,onoff_mic,true_val_mic] = fp_get_signrank_results_megmeg(mic_true,mic_shuf,alpha);
    [true_p_mim,onoff_mim,true_val_mim] = fp_get_signrank_results_megmeg(mim_true,mim_shuf,alpha);
    toc
    
elseif strcmp(testmethod,'t')
    threshold = prctile(reshape(sum(mic_shuf,1),1,[]),99.9);
    true_val_mic = squeeze(sum(mic_true,1));
    onoff_mic = true_val_mic > threshold;
    true_p_mic = [];
    
else
    error('Testmethod unknown.')
end


fprintf('Finding clusters...\n')
kron_conn = fp_get_kron_conn_megmeg(nfreq);
[true_clu_mic, true_total_mic] = fp_get_cluster_components_megmeg(onoff_mic,kron_conn);
[true_clu_mim, true_total_mim] = fp_get_cluster_components_megmeg(onoff_mim,kron_conn);
%findclusters is not appropriate for this application

%% shuffled

for iit = 1:nit
    
    fprintf('Working on iteration %d \n',iit)
    clear onoff cCoh csCoh
    
    %select current shuffled coh
    cCoh_mic = squeeze(mic_shuf(:,iit,:,:,:));
    cCoh_mim = squeeze(mim_shuf(:,iit,:,:,:));
    
    fprintf('Testing...\n')
    tic
    if strcmp(testmethod,'s')
        
        if iit == 1
            csCoh_mic = mic_shuf(:,2:end,:,:,:);
            csCoh_mim = mim_shuf(:,2:end,:,:,:);
        elseif iit == nit
            csCoh_mic = mic_shuf(:,1:end-1,:,:,:);
            csCoh_mim = mim_shuf(:,1:end-1,:,:,:);
        else
            csCoh_mic = mic_shuf(:,[1:iit-1 iit+1:end],:,:,:);
            csCoh_mim = mim_shuf(:,[1:iit-1 iit+1:end],:,:,:);
        end
        
        [~,onoff_mic,shuf_val_mic(iit,:,:,:)] = fp_get_signrank_results_megmeg(cCoh_mic,csCoh_mic,alpha);
        [~,onoff_mim,shuf_val_mim(iit,:,:,:)] = fp_get_signrank_results_megmeg(cCoh_mim,csCoh_mim,alpha);
        
    elseif strcmp(testmethod,'t')
        onoff_mic = squeeze(sum(cCoh_mic,1)) > threshold;
        shuf_val_mic(iit,:,:,:) = squeeze(sum(cCoh_mic,1));
    end
    toc
    
    fprintf('Finding clusters...\n')
    [shuf_clu_mic(iit,:,:,:), shuf_total_mic(iit)] = fp_get_cluster_components_megmeg(onoff_mic,kron_conn);
    [shuf_clu_mim(iit,:,:,:), shuf_total_mim(iit)] = fp_get_cluster_components_megmeg(onoff_mim,kron_conn);
        
    
end


%%
p_mic = fp_get_cluster_p_megmeg(true_total_mic, shuf_total_mic, true_val_mic, shuf_val_mic, true_clu_mic, shuf_clu_mic, fwf);
p_mim = fp_get_cluster_p_megmeg(true_total_mim, shuf_total_mim, true_val_mim, shuf_val_mim, true_clu_mim, shuf_clu_mim, fwf);

%%
outname = sprintf('./mim90_pval');

save(outname,'p_mic','true_total_mic','true_clu_mic','true_p_mic',...
    'true_val_mic','p_mim','true_total_mim','true_clu_mim','true_p_mim',...
    'true_val_mim','-v7.3')
