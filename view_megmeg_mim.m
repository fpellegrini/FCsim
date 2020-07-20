load('ROI_new_j.mat')
load('labels_atlas.mat')

imagesc(squeeze(sum(true_clu,3)))
xticks = linspace(1,28, 28);
% xticklabels = (ROI_new.label)';
yticks = linspace(1,28, 28);
yticklabels = xticklabels;
set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick', yticks, 'YTickLabel',yticklabels)
xtickangle(45)

figure
bar(squeeze(sum(sum(true_clu,1),2)))
xticklabels = 0:5:90;
xticks = linspace(0,45, numel(xticklabels));
set(gca,'XTick', xticks,'XTickLabel',xticklabels)
xlabel('frequency (Hz)')


%%

u = sum(sum(true_clu,2),3);
figure
bar(u)
xticks = linspace(1,28, 28);
set(gca,'XTick', xticks,'XTickLabel',xticklabels)
xtickangle(45)

%%
alpha = [4 6]; 
beta = [7 15];
gamma = [16:45];

imagesc(squeeze(sum(true_clu(:,:,gamma),3)))
xticks = linspace(1,28, 28);
% xticklabels = (ROI_new.label)';
yticks = linspace(1,28, 28);
yticklabels = xticklabels;
set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick', yticks, 'YTickLabel',yticklabels)
xtickangle(45)