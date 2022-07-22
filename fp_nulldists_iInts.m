
for iit = 1:50
    
    cc = triu(abs(randn(68,68)));
    cc = cc+cc';
    a = randperm(68,2);
    iroi_seed= a(1);
    iroi_tar = a(2);
    cc(a(1),a(2)) = 10;
    cc(a(2),a(1)) = 10;
    [mrr(iit), pr(iit),hk(iit),em1(iit),em2(iit),em3(iit)] = fp_mrr_hk(cc,iroi_seed,iroi_tar,1);
    
end


%%

subplot(2,3,1)
[h, u] = fp_raincloud_plot(mrr, [0 0 1], 1,0.2, 'ks');
view([-90 -90]);
set(gca, 'Xdir', 'reverse');
set(gca, 'XLim', [0 1]);
title('Perfect skill MRR')

subplot(2,3,2)
[h, u] = fp_raincloud_plot(pr, [0 0 1], 1,0.2, 'ks');
view([-90 -90]);
set(gca, 'Xdir', 'reverse');
set(gca, 'XLim', [0 1]);
title('Perfect skill PR')

subplot(2,3,3)
[h, u] = fp_raincloud_plot(hk, [0 0 1], 1,0.2, 'ks');
view([-90 -90]);
set(gca, 'Xdir', 'reverse');
set(gca, 'XLim', [0 1]);
title('Perfect skill HK')

subplot(2,3,4)
[h, u] = fp_raincloud_plot(em1, [0 0 1], 1,0.2, 'ks');
view([-90 -90]);
set(gca, 'Xdir', 'reverse');
set(gca, 'XLim', [0 1]);
title('Perfect skill (1-EM1)')

subplot(2,3,5)
[h, u] = fp_raincloud_plot(em2, [0 0 1], 1,0.2, 'ks');
view([-90 -90]);
set(gca, 'Xdir', 'reverse');
set(gca, 'XLim', [0 1]);
title('Perfect skill (1-EM2)')

subplot(2,3,6)
[h, u] = fp_raincloud_plot(em3, [0 0 1], 1,0.2, 'ks');
view([-90 -90]);
set(gca, 'Xdir', 'reverse');
set(gca, 'XLim', [0 1]);
title('Perfect skill (1-EM3)')

%%

load('processed_bs_wzb_90_2000/bs_results.mat')

iatl = 3; % DK atlas
neighbor_thresh = 10; % 10mm vicinity defines neighborhood between regions

[conndist_full, nROI, ROIadj, DconnROI] = get_ROI_dist_full(cortex, iatl, neighbor_thresh);


for iit = 1:50
    
    cc = abs(randn(68,68));
    a = randperm(68,2);
    iroi_seed= a(1);
    iroi_tar = a(2);
    
    
    true_conn = zeros(nROI);
    for iconn = 1:length(iroi_seed)
        true_conn(iroi_seed(iconn), iroi_tar(iconn)) = 1;
    end
    true_conn = true_conn + true_conn';
    
    d1 = true_conn(:);
    d1 = d1 ./ sum(d1(:));
    
    d2 = cc(:);
    d2 = d2 - min(d2);
    
    d2_ = d2;
    d2_ = d2_ ./ sum(d2_(:));
    em1(iit) = 1 - emd_hat_gd_metric_mex(d1, d2_, conndist_full);
    
    d2_ = d2;
    d2_(d2_ < prctile(d2_, 95)) = 0;
    d2_ = d2_ ./ sum(d2_(:));
    em2(iit) = 1 - emd_hat_gd_metric_mex(d1, d2_, conndist_full);
    
    d2_ = d2;
    [so in] = sort(d2_, 'descend');
    d2_(d2_ < so(1)) = 0;
    d2_ = d2_ ./ sum(d2_(:));
    em3(iit) = 1 - emd_hat_gd_metric_mex(d1, d2_, conndist_full);
    
end










