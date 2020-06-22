load('ROI_new_j.mat')

imagesc(squeeze(sum(true_clu,3)))
xticks = linspace(1,28, 28);
xticklabels = (ROI_new.label)';
yticks = linspace(1,28, 28);
yticklabels = (ROI_new.label)';
set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick', yticks, 'YTickLabel',yticklabels)
xtickangle(45)

figure
bar(squeeze(sum(sum(true_clu,1),2)))
xticklabels = 0:5:90;
xticks = linspace(0,45, numel(xticklabels));
set(gca,'XTick', xticks,'XTickLabel',xticklabels)
xlabel('frequency (Hz)')