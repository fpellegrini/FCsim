figure; 

subplot(1,2,1)
a1 = zeros(size(squeeze(true_val(:,1,:))));
a1 = a1-50; 
b = squeeze(true_val(:,1,:));
a1(squeeze(p_masked_h(:,1,:)))= b(squeeze(p_masked_h(:,1,:)));
a1(squeeze(p_masked_l(:,1,:)))= b(squeeze(p_masked_l(:,1,:)));

imagesc(a1)
title('Right fdr-masked test statistic')
xlabel('Frequency in Hz') 
ylabel('Voxels')
xticklabels = 0:5:90;
xticks = linspace(0,45, numel(xticklabels));
set(gca,'XTick', xticks,'XTickLabel',xticklabels)
colorbar

subplot(1,2,2)
a2 = zeros(size(squeeze(true_val(:,2,:))));
a2 = a2-50; 
b = squeeze(true_val(:,2,:));
a2(squeeze(p_masked_h(:,2,:)))= b(squeeze(p_masked_h(:,2,:)));
a2(squeeze(p_masked_l(:,2,:)))= b(squeeze(p_masked_l(:,2,:)));

imagesc(a2)
title('Left fdr-masked test statistic') 
xlabel('Frequency in Hz') 
ylabel('Voxels')
xticklabels = 0:5:90;
xticks = linspace(0,45, numel(xticklabels));
set(gca,'XTick', xticks,'XTickLabel',xticklabels)
colorbar

%%
o = [0 130];
ch = squeeze(sum(p_masked_h,1));
cl = squeeze(sum(p_masked_l,1));

figure
subplot(2,2,1)
bar([0:2:90],squeeze(ch(1,:)))
title('positive fdr significant gc, right lfp')
ylim(o)
xlabel('frequency in Hz')
subplot(2,2,2)
bar([0:2:90],squeeze(ch(2,:)))
title('positive fdr significant gc, left lfp ')
ylim(o)
xlabel('frequency in Hz')
subplot(2,2,3)
bar([0:2:90],squeeze(cl(1,:)))
title('negative fdr significant gc, right lfp')
ylim(o)
xlabel('frequency in Hz')
subplot(2,2,4)
bar([0:2:90],squeeze(cl(2,:)))
title('negative fdr significant gc, left lfp')
ylim(o)
xlabel('frequency in Hz')

%%

dh = squeeze(sum(p_masked_h,3));
dl = squeeze(sum(p_masked_l,3)); 

outname = sprintf('cluster_gc_01_%s%s_%d_lcmv_welch.nii',sides(iside), vals(ival),iclus);
fp_data2nii(b,pos,[],outname)





