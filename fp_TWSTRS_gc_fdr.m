function fp_TWSTRS_gc_fdr

DIROUT = './';
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[~, voxID] = fp_find_commonvox;

scores = db_twstrs();
alpha = 0.005;
nit = 1000;
%%
load(sprintf('%sDIFFGC_lcmv.mat',DIROUT));
[nsubs,nvox,nside,nfreq] = size(DIFFGC);
%%
o = 1;
for id = 1:nsubs
    if ~isnan(scores{1,id}{1,2})
        fprintf('Working on subject %s \n',patientID{id})
        
        t_score(o) = scores{1,id}{1,2};
        diffgc_t(o,:,:,:) = DIFFGC(id,:,:,:);
        o = o+1;
        
    end
end
clear o
%% true

for iside = 1:2 
    for ifreq = 1:nfreq 
        for ivox = 1:nvox 
            [rho(ivox,ifreq, iside), pval(ivox, ifreq, iside)] ...
                = corr(t_score', diffgc_t(:,ivox,iside,ifreq), 'Type','Pearson');
        end 
    end 
end

[p, mask] = fdr(pval,0.05); 


