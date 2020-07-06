function fp_gc_fdr(DIROUT, alpha,fwf,j)

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
    [~, voxID] = fp_find_commonvox;
elseif j == 1
    j_s = 'j';
    patientID = {'04'; '07'; '09'; '10';'11';'20';'22';'25'};
    [~, voxID] = fp_find_commonvox;
    voxID([3,7,8])=[];
else
    error('Wrong input j')
end


load(sprintf('%sDIFFGC_lcmv',DIROUT));
[nsubs,nvox,nside,nfreq] = size(DIFFGC);
nit = 500;

%% true
fprintf('Testing...\n')
tic
[true_p,onoff,true_val,true_effectdir] = fp_get_signrank_results_gc(DIFFGC,alpha);
toc

[p_fdr_true, p_masked_true] = fdr( true_p, 0.05);

bh =zeros(size(true_val)); 
bh(true_effectdir>0) = true_val(true_effectdir>0);
bl =zeros(size(true_val)); 
bl(true_effectdir<0) = true_val(true_effectdir<0);
%% shuffled

for iit = 1:nit
    
    fprintf('Working on iteration %d \n',iit)
    clear onoff c_DIFFGC effectdir p_onoff n_onoff
    
    %shuffle signs of DIFFGC
    c_DIFFGC = DIFFGC .* (sign(randn(size(DIFFGC))));
    
    fprintf('Testing...\n')
    [~,~,shuf_val, shuf_effectdir] = fp_get_signrank_results_gc(c_DIFFGC,alpha);
    
    al = zeros(size(shuf_val)); 
    ah = al; 
    al(shuf_effectdir<0) = shuf_val(shuf_effectdir<0);
    ah(shuf_effectdir>0) = shuf_val(shuf_effectdir>0);        
    
    uh(iit,:,:,:) = ah>bh;
    ul(iit,:,:,:) = al<bl;
    
end

ph = squeeze(sum(uh,1))./nit;
pl = squeeze(sum(ul,1))./nit;

[p_fdr_l, p_masked_l] = fdr( pl, 0.025);
[p_fdr_h, p_masked_h] = fdr( ph, 0.025);


outname=[DIROUT 'gc_fdr.mat'];
save(outname,'p_fdr_true','p_masked_true','p_fdr_l','p_masked_l','p_fdr_h','p_masked_h','true_val','true_effectdir')
