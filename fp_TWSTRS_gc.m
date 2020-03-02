function fp_TWSTRS_gc

DIROUT = './';
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[~, voxID] = fp_find_commonvox;

scores = db_twstrs();
alpha = 0.05;
nit = 1000;
%%
load(sprintf('%sDIFFGC.mat',DIROUT));
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
            [true_rho(ivox,ifreq, iside), true_pval(ivox, ifreq, iside)] ...
                = corr(t_score', diffgc_t(:,ivox,iside,ifreq), 'Type','Spearman');
        end 
    end 
end

kron_conn = fp_get_kron_conn(patientID{id}, voxID{id}, nfreq); %same for every id
for iside = 1:2 
    clear onoff onoff1 
    onoff = true_pval(:,:,iside) < alpha;

    %pos
    onoff1 = onoff;
    onoff1(true_rho(:,:,iside)<0)=0;
    [true_clu_pos(:,:,iside), true_total_pos(iside)] = fp_get_cluster_components(onoff1,kron_conn);

    %neg
    onoff1 = onoff;
    onoff1(true_rho(:,:,iside)>0)=0;
    [true_clu_neg(:,:,iside), true_total_neg(iside)] = fp_get_cluster_components(onoff1,kron_conn);
end

%% shuffled

for iit = 1:nit
    fprintf('Working on iteration %d \n',iit);
    clear pval_shuf onoff onoff1 
    
    diffgc_s = diffgc_t .* (sign(randn(size(diffgc_t)))); 
    
    for iside = 1:2
        clear onoff pval_shuf 
        for ifreq = 1:nfreq
            for ivox = 1:nvox
                [rho_shuf(iit,ivox,ifreq,iside), pval_shuf(ivox, ifreq)] = corr(t_score', diffgc_s(:,ivox,iside,ifreq), 'Type','Spearman');
            end
        end
        onoff = pval_shuf < alpha;
        
        %pos
        onoff1 = onoff;
        onoff1(squeeze(rho_shuf(iit,:,:,iside))<0)=0;
        [shuf_clu(iit,:,:,1, iside), shuf_total(iit,1,iside)] = fp_get_cluster_components(onoff1,kron_conn);
        
        %neg
        onoff1 = onoff;
        onoff1(squeeze(rho_shuf(iit,:,:,iside))>0)=0;
        [shuf_clu(iit,:,:,2, iside), shuf_total(iit,2,iside)] = fp_get_cluster_components(onoff1,kron_conn);
    end
end
    


%%

for iside = 1:2
    p_pos{iside} = fp_get_cluster_p_gc_new(true_total_pos(iside), shuf_total(:,:,iside),...
        true_rho(:,:,iside), rho_shuf(:,:,:,iside), true_clu_pos(:,:,iside), shuf_clu(:,:,:,:,iside), 0);
    p_neg{iside} = fp_get_cluster_p_gc_new(true_total_neg(iside), shuf_total(:,:,iside),...
        true_rho(:,:,iside), rho_shuf(:,:,:,iside), true_clu_neg(:,:,iside), shuf_clu(:,:,:,:,iside), 0);
end

outname = './TWSTRS_gc.mat';
save(outname,'p_pos','p_neg','true_total_pos','true_total_neg','true_clu_pos','true_clu_neg','true_pval','true_rho','-v7.3')
