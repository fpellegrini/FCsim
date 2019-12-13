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


nsubs = numel(patientID);
minnbchan = 2;

    load(sprintf('%sDIFFGC',DIROUT));


[nsubs,nvox,nside,nfreq] = size(DIFFGC);


%% true
fprintf('Testing...\n')
tic
[true_p,onoff,true_val] = fp_get_signrank_results_gc(DIFFGC,alpha);
toc

fprintf('Finding clusters...\n')
kron_conn = fp_get_kron_conn_gc(nfreq);

for iside = 1:nside
    clear c_onoff
    c_onoff = squeeze(onoff(:,:,iside,:));
    [true_clu(:,:,iside), true_total(iside)] = fp_get_cluster_components_megmeg(c_onoff,kron_conn);
end

%% shuffled

for iit = 1:nit
    
    fprintf('Working on iteration %d \n',iit)
    clear onoff c_DIFFGC
    
    %shuffle signs of DIFFGC
    c_DIFFGC = DIFFGC .* reshape(sign(DIFFGC(randperm(numel(DIFFGC)))),size(DIFFGC));
    
    fprintf('Testing...\n')
    [~,onoff,shuf_val(iit)] = fp_get_signrank_results_gc(c_DIFFGC,alpha);
    
    for iside=1:nside
        clear c_onoff
        c_onoff = squeeze(onoff(:,:,iside,:));
        [shuf_clu(iit,:,:,iside), shuf_total(iit,iside)] = fp_get_cluster_components_gc(c_onoff,kron_conn);
    end
    
end
%%
%%%%%%% test this with data, not ready yet! 
for iside=1:nside
    p(iside) = fp_get_cluster_p_megmeg(true_total, shuf_total, true_val, shuf_val, true_clu, shuf_clu, fwf);
end

%%
outname = sprintf('%sp_gc_%s_%s_%s',DIROUT, alpha_s, fwf_s, j_s);

save(outname,'p','true_total','true_clu','true_p','true_val','-v7.3')
