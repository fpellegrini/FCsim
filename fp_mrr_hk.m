function [mrr, mrrs,hk] = fp_mrr_hk(cc,iroi_seed,iroi_tar)

if norm(iroi_seed(:)- iroi_tar(:),'fro')==0
    mrr = 0;
    mrrs=0;
    hk = 0; 
else

    nints = numel(iroi_seed);
    nroi = size(cc,1);
    inds = find(triu(ones(nroi),1));

    label= zeros(nroi,nroi);
    for iint = 1:numel(iroi_seed)
        label(iroi_seed(iint),iroi_tar(iint)) = 1;
        label(iroi_tar(iint),iroi_seed(iint)) = 1;
    end
    lab = label(:);
    tr = find(lab(inds)==1);
    k=10;

    [ss, idx] = sort(cc(inds),'descend');

    for it = 1:numel(tr)
        r1(it) = find(idx==tr(it));
    end

    mrr = sum(1./r1)/nints;
    mrrs = 1-sum(r1/k);
    hk = sum(r1<=k)/nints;
end