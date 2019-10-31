function fp_cluster_proto(DIROUT,clusterfun,abs_imag,testmethod,alpha,fwf, j)
%Group statistics.
%Clustering group statistics across space and frequencies.
%
%Input:
%
%clusterfun = either 'c' for components fun or 'f' for findclusters fun 
%
%abs_imag = either 'abs' or 'imag'
%
%testmethod = either 's' for signrank testing or 't' for simple
%thresholding
%
%fwf: testing method (less conservative): test first cluster against
%first, second with second etc. By default, clusters are always compared
%against the largest shuffled cluster. 
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
    [~, voxID] = fp_find_commonvox;
elseif j_s == 1
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
%%

for id = 1:nsubs    
    fprintf('Working on subject %s',patientID{id})
    
    %get neighbouring nodes and node positions
    clear flip_id noEq
    
    %get flip id and symmetric head
    [flip_id, noEq] = fp_get_flip_id(patientID{id},voxID{id});
    
    for ichunk = 1:nchunk
 
        inname = sprintf('Coherences_e_Patient%s_chunk%d.mat',patientID{id},ichunk);                
        COH1(id,ichunk,:,:,:) = fp_get_coh(inname, noEq, vox_ind,flip_id, abs_imag);
        
    end   
end

[nsub, nchunk,niit,nfreq,ns] = size(COH1);
COH = reshape(COH1,[nsub,nchunk*niit,nfreq,ns]);
nit = size(COH,2)-1; %update nit to number of *shuffled* it

%true and shuffled coherences
tCoh = squeeze(COH(:,1,:,:));
sCoh = COH(:,2:end,:,:);

%% true 
if strcmp(testmethod,'s')    
    [true_p,onoff,true_val] = fp_get_signrank_results(tCoh,sCoh,alpha);
    
elseif strcmp(testmethod,'t')
    threshold = prctile(reshape(sCoh,1,[]),99.9);
    onoff = squeeze(sum(tCoh,1)) > threshold;
    true_val = squeeze(sum(tCoh,1));
    true_p = [];
    
else
    error('Testmethod unknown.')
end

if strcmp(clusterfun, 'c')
    kron_conn = fp_get_kron_conn(patientID{id}, voxID{id}, nfreq); %same for every id
    [true_clu, true_total] = fp_get_cluster_components(onoff,kron_conn);
    
elseif strcmp(clusterfun, 'f')
    conn = fp_find_neighbours(patientID{id});
    match_conn = conn(voxID{id},voxID{id}); %same for evey id
    [true_clu, true_total] = fp_get_cluster_findclusters(onoff, match_conn,minnbchan);
    
else
    error('Cluster method unknown.')
end

%% shuffled

for iit = 1:nit 
    clear onoff cCoh csCoh
    
    %select current shuffled coh
    cCoh = squeeze(sCoh(:,iit,:,:)); 
    
    if strcmp(testmethod,'s')
        
        if iit == 1
            csCoh = sCoh(:,2:end,:,:);
        elseif iit == nit
            csCoh = sCoh(:,1:end-1,:,:);
        else
            csCoh = sCoh(:,[1:iit-1 iit+1:end],:,:);
        end
        
        [~, onoff, shuf_val(iit,:,:)] = fp_get_signrank_results(cCoh,csCoh,alpha);
        
    elseif strcmp(testmethod,'t')
        onoff = squeeze(sum(cCoh,1)) > threshold;
        shuf_val(iit,:,:) = squeeze(sum(cCoh,1));
    end
    
    if strcmp(clusterfun, 'c')
        [shuf_clu(iit,:,:), shuf_total(iit)] = fp_get_cluster_components(onoff,kron_conn);
        
    elseif strcmp(clusterfun, 'f')
        [shuf_clu(iit,:,:), shuf_total(iit)] = fp_get_cluster_findclusters(onoff, match_conn,minnbchan);
        
    else
    end
    
end
%%
p = fp_get_cluster_p(true_total, shuf_total, true_val, shuf_val, true_clu, shuf_clu, fwf);

%%
outname = sprintf('%sp_cluster_g_%s_freq_%s_%s_%s_%s_%s',DIROUT,clusterfun, abs_imag,...
    testmethod, alpha_s, fwf_s, j_s);

save(outname,'p','true_total','true_clu','true_p','true_val','-v7.3')
