function fp_cluster_megmeg(DIROUT, abs_imag,testmethod,alpha,fwf,j)

%Group statistics on megmeg data.
%Clustering group statistics across space and frequencies.
%
%Input:
%
%abs_imag = either 'abs' or 'imag'
%
%testmethod = either 's' for signrank testing or 't' for simple
%thresholding
%
%fwf: testing method (less conservative): test first cluster against
%
%first, second with second etc. By default, clusters are always compared
%against the largest shuffled cluster.
%
%j=1 : only Julian's subjects


% fp_addpath_sabzi

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
    [~, voxID] = fp_find_commonvox;
elseif j == 1
    j_s = 'j';
    patientID = {'04'; '07'; '09'; '10';'11';'20';'22';'25'};
    [~, voxID] = fp_find_commonvox;
    voxID([3,7,8])=[];
else
    error('Wrong input j')
end


nsubs = numel(patientID);
nchunk = 50;
minnbchan = 2;


for id = 1:5%nsubs
    fprintf('Working on subject %s \n',patientID{id})
    
    clear COH true_coh r_coh rt_coh
    load(sprintf('%sroi_coh_sub%s',DIROUT,patientID{id}));
    
    %absolute value
    if strcmp(abs_imag,'abs')
        r_coh= abs(COH);
        rt_coh = abs(true_coh);
    elseif strcmp(abs_imag,'imag')
        r_coh = abs(imag(COH));
        rt_coh = abs(imag(true_coh));
    else
        error('Method unknown!')
    end
    
    sCoh(id,:,:,:,:) = r_coh;
    tCoh(id,:,:,:) = rt_coh;
    
end

[nsubs,nit,nroi,~,nfreq] = size(sCoh);

%% true
fprintf('Testing...\n')
if strcmp(testmethod,'s')
    %     [true_p,onoff,true_val] = fp_get_signrank_results(tCoh,sCoh,alpha);
    
elseif strcmp(testmethod,'t')
    threshold = prctile(reshape(sum(sCoh,1),1,[]),99.9);
    true_val = squeeze(sum(tCoh,1));
    onoff = true_val > threshold;
    true_p = [];
    
else
    error('Testmethod unknown.')
end


fprintf('Finding clusters...\n')
kron_conn = fp_get_kron_conn_megmeg(nfreq);
[true_clu, true_total] = fp_get_cluster_components_megmeg(onoff,kron_conn);
%findclusters is not appropriate for this application

%% shuffled

for iit = 1:nit
    
    fprintf('Working on iteration %d \n',iit)
    clear onoff cCoh csCoh
    
    %select current shuffled coh
    cCoh = squeeze(sCoh(:,iit,:,:,:));
    
    fprintf('Testing...\n')
    tic
    if strcmp(testmethod,'s')
        
        %         if iit == 1
        %             csCoh = sCoh(:,2:end,:,:);
        %         elseif iit == nit
        %             csCoh = sCoh(:,1:end-1,:,:);
        %         else
        %             csCoh = sCoh(:,[1:iit-1 iit+1:end],:,:);
        %         end
        %
        %         [~, onoff, shuf_val(iit,:,:)] = fp_get_signrank_results(cCoh,csCoh,alpha);
        %
    elseif strcmp(testmethod,'t')
        onoff = squeeze(sum(cCoh,1)) > threshold;
        shuf_val(iit,:,:,:) = squeeze(sum(cCoh,1));
    end
    toc
    
    fprintf('Finding clusters...\n')
    [shuf_clu(iit,:,:,:), shuf_total(iit)] = fp_get_cluster_components_megmeg(onoff,kron_conn);
    
    
end
%%
p = fp_get_cluster_p_megmeg(true_total, shuf_total, true_val, shuf_val, true_clu, shuf_clu, fwf);

%%
outname = sprintf('%sp_cluster_g_freq_%s_%s_%s_%s_%s_%s',DIROUT, abs_imag,...
    testmethod, alpha_s, fwf_s, j_s, filtertype);

save(outname,'p','true_total','true_clu','true_p','true_val','-v7.3')
