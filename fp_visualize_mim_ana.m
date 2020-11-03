
clear all
close all

load('./mimsim_ana_snr02_30nit.mat')

figure
bar(median_auc)
title('median auc')
xlabel('pipelines')
ylabel('auc')
legend('mic','mim')
xticks = 1:9; 
xticklabels = {'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline'};
set(gca,'XTick', xticks,'XTickLabel',xticklabels)
figure
bar(-log10(p_auc))
title('-log10(p auc)')
xlabel('pipelines')
ylabel('-log10(p auc)')
legend('mic','mim')
xticks = 1:9; 
xticklabels = {'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline'};
set(gca,'XTick', xticks,'XTickLabel',xticklabels)

figure
subplot(2,1,1)
boxplot(squeeze(perf_mim_corr'),'Labels',{'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline'})
grid on
title('MIM correlation with ground truth')
xlabel('Pipelines')
ylabel('correlation coefficient')

subplot(2,1,2)
boxplot(squeeze(perf_mic_corr'),'Labels',{'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline'})
grid on
title('MIC correlation with ground truth')
xlabel('Pipelines')
ylabel('correlation coefficient')