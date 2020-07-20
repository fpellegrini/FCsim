function fp_gc_fdr(DIROUT)


patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};


load(sprintf('./DIFFGC_lcmv'));
[nsubs,nvox,nside,nfreq] = size(DIFFGC);
nit = 500;

%% true
fprintf('Testing...\n')
tic
[true_p,onoff,true_val,true_effectdir] = fp_get_signrank_results_gc(DIFFGC,0.005);
toc

[p, mask] = fdr( true_p, 0.05);
