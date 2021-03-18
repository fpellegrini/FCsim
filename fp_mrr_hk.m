function [mrr, pr,hk] = fp_mrr_hk(cc,iroi_seed,iroi_tar,half_m)

if norm(iroi_seed(:)- iroi_tar(:),'fro')==0
    mrr = 0;
    pr=0;
    hk = 0; 
else

    nints = numel(iroi_seed);
    nroi = size(cc,1);
    
    if half_m 
        inds = find(triu(ones(nroi),1));
    else
        inds = find(ones(nroi));
    end

    label= zeros(nroi,nroi);
    for iint = 1:numel(iroi_seed)
        label(iroi_seed(iint),iroi_tar(iint)) = 1;
        if half_m %when only positive DIFFGC is submitted as cc
            label(iroi_tar(iint),iroi_seed(iint)) = 1;
        end
    end
    lab = label(:);
    tr = find(lab(inds)==1);
    k=10;
    for ii = 1:nints
       perfectSkill(ii) = 1/ii;
       noSkill(ii) = 1/(numel(lab(inds))-ii+1);
    end
    perfectSkill = sum(perfectSkill)/nints;
    noSkill = sum(noSkill)/nints;

    [ss, idx] = sort(cc(inds),'descend');

    for it = 1:numel(tr)
        r1(it) = find(idx==tr(it));
    end
    
    mrr = sum(1./r1)/nints;
    mrr = (mrr-noSkill)/(perfectSkill-noSkill);
    
    hk = sum(r1<=k)/nints;
    
    if any(r1>nints)
        [~,~,~,pr]=perfcurve(lab(inds),cc(inds),1,'XCrit','reca','YCrit','prec');
        
    else %perfect skill
        pr=1;
    end
    
    
    
end