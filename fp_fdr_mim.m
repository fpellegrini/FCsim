function fp_fdr_mim(DIROUT,alpha,fwf,j)

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


fp_addpath_sabzi

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
[true_p_mic,onoff_mic,true_val_mic] = fp_get_signrank_results_megmeg(mic_true,mic_shuf,alpha);
[true_p_mim,onoff_mim,true_val_mim] = fp_get_signrank_results_megmeg(mim_true,mim_shuf,alpha);

[p_fdr_true_mic, p_masked_true_mic] = fdr( true_p_mic, 0.05);
[p_fdr_true_mim, p_masked_true_mim] = fdr( true_p_mim, 0.05);
%% shuffled

for iit = 1:nit
   
    
    
    fprintf('Working on iteration %d \n',iit)
    clear onoff cCoh csCoh
    
    %select current shuffled coh
    cCoh_mic = squeeze(mic_shuf(:,iit,:,:,:));
    cCoh_mim = squeeze(mim_shuf(:,iit,:,:,:));
    
    fprintf('Testing...\n')
    
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
    
    [shuf_p_mic(iit,:,:,:),~,shuf_val_mic(iit,:,:,:)] = fp_get_signrank_results_megmeg(cCoh_mic,csCoh_mic,alpha);
    [shuf_p_mim(iit,:,:,:),~,shuf_val_mim(iit,:,:,:)] = fp_get_signrank_results_megmeg(cCoh_mim,csCoh_mim,alpha);
    
end


for iit = 1:nit 
     u_mic(iit,:,:,:) = squeeze(shuf_val_mic(iit,:,:,:))>true_val_mic;
     u_mim(iit,:,:,:) = squeeze(shuf_val_mim(iit,:,:,:))>true_val_mic;
end

p_mic = squeeze(sum(u_mic,1))./nit;
p_mim = squeeze(sum(u_mim,1))./nit;

[p_fdr_mic, p_masked_mic] = fdr( p_mic, 0.05);
[p_fdr_mim, p_masked_mim] = fdr( p_mim, 0.05);


outname=[DIROUT 'megmeg_fdr.mat'];
save(outname,'p_fdr_true_mic','p_fdr_true_mim','p_masked_true_mic','p_masked_true_mim',...
    'p_fdr_mic','p_fdr_mim','p_masked_mic','p_masked_mim','true_val_mim','true_val_mic')

