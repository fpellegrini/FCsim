function fp_cluster_gc(DIROUT, alpha,fwf,j)

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


load(sprintf('%sDIFFGC',DIROUT));
[nsubs,nvox,nside,nfreq] = size(DIFFGC);
nit = 1000;

%% true
fprintf('Testing...\n')
tic
[true_p,onoff,true_val,effectdir] = fp_get_signrank_results_gc(DIFFGC,alpha);
toc

fprintf('Finding clusters...\n')
kron_conn = fp_get_kron_conn_gc(nfreq, voxID);

%pos 
p_onoff= onoff;
p_onoff(effectdir<0)=0;

for iside = 1:nside
    clear c_onoff
    c_onoff = squeeze(p_onoff(:,iside,:));    
    [true_clu(:,:,iside,1), true_total(iside,1)] = fp_get_cluster_components_gc(c_onoff,kron_conn);
end

%neg 
n_onoff= onoff;
n_onoff(effectdir>0)=0;

for iside = 1:nside
    clear c_onoff
    c_onoff = squeeze(n_onoff(:,iside,:));
    
    [true_clu(:,:,iside,2), true_total(iside,2)] = fp_get_cluster_components_gc(c_onoff,kron_conn);
end

%% shuffled

for iit = 1:nit
    
    fprintf('Working on iteration %d \n',iit)
    clear onoff c_DIFFGC effectdir p_onoff n_onoff
    
    %shuffle signs of DIFFGC
    c_DIFFGC = DIFFGC .* (sign(randn(size(DIFFGC))));
    
    fprintf('Testing...\n')
    [~,onoff,shuf_val(iit,:,:,:), effectdir] = fp_get_signrank_results_gc(c_DIFFGC,alpha);
    
    %pos
    p_onoff= onoff;
    p_onoff(effectdir<0)=0;
    
    for iside = 1:nside
        clear c_onoff
        c_onoff = squeeze(p_onoff(:,iside,:));
        
        [shuf_clu(iit,:,:,iside,1), shuf_total(iit,iside,1)] = fp_get_cluster_components_gc(c_onoff,kron_conn);
    end
    
    %neg
    n_onoff= onoff;
    n_onoff(effectdir>0)=0;
    
    for iside = 1:nside
        clear c_onoff
        c_onoff = squeeze(n_onoff(:,iside,:));
        
        [shuf_clu(iit,:,:,iside,2), shuf_total(iit,iside,2)] = fp_get_cluster_components_gc(c_onoff,kron_conn);
    end   
    
end
%%
for iside=1:nside
    %pos
    p{iside,1} = fp_get_cluster_p_gc_new(true_total(iside,1), shuf_total(:,iside,:), squeeze(true_val(:,iside,:)),...
        squeeze(shuf_val(:,:,iside,:)), true_clu(:,:,iside,1), shuf_clu(:,:,:,iside,:), fwf);
    %neg
    p{iside,2} = fp_get_cluster_p_gc_new(true_total(iside,2), shuf_total(:,iside,:), squeeze(true_val(:,iside,:)),...
        squeeze(shuf_val(:,:,iside,:)), true_clu(:,:,iside,2), shuf_clu(:,:,:,iside,:), fwf);
end

%%
outname = sprintf('%sp_gc_%s_%s_%s',DIROUT, alpha_s, fwf_s, j_s);

save(outname,'p','true_total','true_clu','true_p','true_val','-v7.3')
