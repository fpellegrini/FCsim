function [mrr, pr,hk,em1,em2,em3] = fp_mrr_hk(cc,iroi_seed,iroi_tar,half_m)

if norm(iroi_seed(:)- iroi_tar(:),'fro')==0
    mrr = 0;
    pr = 0;
    hk = 0; 
    em1 = 0;
    em2 = 0;
    em3 = 0;
else

    nints = numel(iroi_seed);
    nroi = size(cc,1);
    
    if half_m 
        ind_ = find(triu(ones(nroi),1));
    else
        ind_ = find(ones(nroi));
    end

    label= zeros(nroi,nroi);
    for iint = 1:numel(iroi_seed)
        label(iroi_seed(iint),iroi_tar(iint)) = 1;
        if half_m %when only positive DIFFGC is submitted as cc
            label(iroi_tar(iint),iroi_seed(iint)) = 1;
        end
    end
    lab = label(:);
    tr = find(lab(ind_)==1);
    
    %% MRR
    for ii = 1:nints
       perfectSkill(ii) = 1/ii;
       noSkill(ii) = 1/(numel(lab(ind_))-ii+1);
    end
    perfectSkill = sum(perfectSkill)/nints;
    noSkill = sum(noSkill)/nints;

    [ss, idx] = sort(cc(ind_),'descend');

    for it = 1:numel(tr)
        r1(it) = find(idx==tr(it));
    end
    
    mrr = sum(1./r1)/nints;
    mrr = (mrr-noSkill)/(perfectSkill-noSkill);
    
    %% hits @ k 
    k=10;
    hk = sum(r1<=k)/nints;
    
    %% area under the precision-recall curve 
    if any(r1>nints)
        [~,~,~,pr]=perfcurve(lab(ind_),cc(ind_),1,'XCrit','reca','YCrit','prec');
        
    else %perfect skill
        pr=1;
    end
    
    %% moving earth distance in 3 variants
    try
        load('processed_bs_wzb_90_2000/bs_results.mat')
        iatl = 3; % DK atlas
        neighbor_thresh = 10; % 10mm vicinity defines neighborhood between regions
%         [conndist, ~,~, ~] = get_ROI_dist2(cortex, iatl, neighbor_thresh);
        [conndist_full, ~, ~, ~] = get_ROI_dist_full(cortex, iatl, neighbor_thresh);

        d1 = label(:); % takes always the full matrix, not only upper half 
        d1 = d1 ./ sum(d1(:));

        d2 = cc(:);
        d2 = d2 - min(d2);

        d2_ = d2;
        d2_ = d2_ ./ nansum(d2_(:));
        em1 = 1 - emd_hat_gd_metric_mex(d1, d2_, conndist_full);

        d2_ = d2;
        d2_(d2_ < prctile(d2_, 95)) = 0;
        d2_ = d2_ ./ nansum(d2_(:));
        em2 = 1 - emd_hat_gd_metric_mex(d1, d2_, conndist_full);

        d2_ = d2;
        [so, ~] = sort(d2_, 'descend');
        d2_(d2_ < so(numel(iroi_seed))) = 0;
        d2_ = d2_ ./ nansum(d2_(:));
        em3 = 1 - emd_hat_gd_metric_mex(d1, d2_, conndist_full);
        
    catch %that not the full simulation stops when compilation of mex files fail 
        
        em1 = 0; 
        em2 = 0; 
        em3 = 0; 
        
    end
    
    
    
end