function fp_mimsim_ana

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim_snr02/';

%default sparamenters
nit = 100;
iInt = 1; 
iReg=1; 
isnr=0.2;
iss = 0.5;
ilag=2;
ihemi=0;
ifilt='l';
%%
for iit=[12 14 18 22 23 8]
    inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d_iter%d'...
        ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,ihemi,iit);
    
    load([DIRIN inname '.mat'])
        
    %ground truth
%     gt.mic = abs(imag(gt.mic));
%     gt.mim = abs(imag(gt.mim)); 
    
    %PERFORMANCE has the dimensions: 2 (mim/mic) x 8 (pipelines)
%     clear PERFORMANCE BASELINE
%     [PERFORMANCE, BASELINE] = fp_get_performance(gt, mic, mim, mean_coh);
% 
%     perf_mim_corr(:,iit) = [squeeze(PERFORMANCE(1,:)) squeeze(BASELINE(1,1))];
%     perf_mic_corr(:,iit) = [squeeze(PERFORMANCE(2,:)) squeeze(BASELINE(1,1))];
%     perf_mc_corr(:,iit) = [squeeze(PERFORMANCE(3,:))];
    
    %AUC
    label= zeros(68,68); 
    label(iroi_seed,iroi_tar) = 1;
    label(iroi_tar,iroi_seed) = 1;
    
    for ipip = 1:5
        cc = sum(mic.fixed{ipip},3);
        [~, ~,~, auc(iit,ipip,1)] = ...
            perfcurve(label(:),cc(:),1);
        clear cc
        cc = sum(mim.fixed{ipip},3);
        [~,~,~, auc(iit,ipip,2)] = ...
            perfcurve(label(:),cc(:),1);
        clear cc
        cc = sum(mean_coh.fixed{ipip},3);
        [~,~,~, auc(iit,ipip,3)] = ...
            perfcurve(label(:),cc(:),1);
        clear cc
        
        corrs(iit,ipip,1) = to_save.fixed{ipip}.corr_voxmic;
        corrs(iit,ipip,2) = to_save.fixed{ipip}.corr_voxmim; 
        corrs(iit,ipip,3) = to_save.fixed{ipip}.corr_voxmeancoh;
        corrs(iit,ipip,4) = corr(to_save.fixed{ipip}.npcs,to_save.nvoxroi');
        
        varex(iit,ipip) = to_save.fixed{ipip}.var_explained;
        
    end
    cc = sum(mic.max,3);
    [~,~,~, auc(iit,6,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.max,3);
    [~,~,~, auc(iit,6,2)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mean_coh.max,3);
    [~,~,~, auc(iit,6,3)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    corrs(iit,6,1) = to_save.max.corr_voxmic;
    corrs(iit,6,2) = to_save.max.corr_voxmim; 
    corrs(iit,6,3) = to_save.max.corr_voxmeancoh;
    corrs(iit,6,4) = corr(to_save.max.npcs',to_save.nvoxroi');
    
    cc = sum(mic.percent,3);
    [~,~,~, auc(iit,7,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.percent,3);
    [~,~,~, auc(iit,7,2)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc 
    cc = sum(mean_coh.percent,3);
    [~,~,~, auc(iit,7,3)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc 
    corrs(iit,7,1) = to_save.percent.corr_voxmic;
    corrs(iit,7,2) = to_save.percent.corr_voxmim; 
    corrs(iit,7,3) = to_save.percent.corr_voxmeancoh;
    corrs(iit,7,4) = corr(to_save.percent.npcs',to_save.nvoxroi');
    
    cc = sum(mic.case2,3);
    [~,~,~, auc(iit,8,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.case2,3);
    [~,~,~, auc(iit,8,2)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    corrs(iit,8,1) = to_save.case2.corr_voxmic;
    corrs(iit,8,2) = to_save.case2.corr_voxmim; 
    
    cc = sum(mic.baseline,3);
    [~,~,~, auc(iit,9,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.baseline,3);
    [~,~,~, auc(iit,9,2)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    corrs(iit,9,1) = to_save.baseline.corr_voxmic;
    corrs(iit,9,2) = to_save.baseline.corr_voxmim; 
    
    cc = sum(mic.max_corrected,3);
    [~,~,~, auc(iit,10,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.max_corrected,3);
    [~,~,~, auc(iit,10,2)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mean_coh.max_corrected,3);
    [~,~,~, auc(iit,10,3)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    corrs(iit,10,1) = to_save.max_corrected.corr_voxmic;
    corrs(iit,10,2) = to_save.max_corrected.corr_voxmim; 
    corrs(iit,10,3) = to_save.max_corrected.corr_voxmeancoh;
    corrs(iit,10,4) = corr(to_save.max_corrected.npcs',to_save.nvoxroi');
    
    cc = sum(mic.percent_corrected,3);
    [~,~,~, auc(iit,11,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.percent_corrected,3);
    [~,~,~, auc(iit,11,2)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc 
    cc = sum(mean_coh.percent_corrected,3);
    [~,~,~, auc(iit,11,3)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc 
    corrs(iit,11,1) = to_save.percent_corrected.corr_voxmic;
    corrs(iit,11,2) = to_save.percent_corrected.corr_voxmim; 
    corrs(iit,11,3) = to_save.percent_corrected.corr_voxmeancoh;
    corrs(iit,11,4) = corr(to_save.percent_corrected.npcs',to_save.nvoxroi');
    
    
    for ipip = 1:5
        cc = sum(mic.fixed_zs0{ipip},3);
        [~, ~,~, auc(iit,11+ipip,1)] = ...
            perfcurve(label(:),cc(:),1);
        clear cc
        cc = sum(mim.fixed_zs0{ipip},3);
        [~,~,~, auc(iit,11+ipip,2)] = ...
            perfcurve(label(:),cc(:),1);
        clear cc
        cc = sum(mean_coh.fixed_zs0{ipip},3);
        [~,~,~, auc(iit,11+ipip,3)] = ...
            perfcurve(label(:),cc(:),1);
        clear cc
        
       
    end
     
    cc = sum(mic.max_zs0,3);
    [~,~,~, auc(iit,17,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.max_zs0,3);
    [~,~,~, auc(iit,17,2)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mean_coh.max_zs0,3);
    [~,~,~, auc(iit,17,3)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    
    cc = sum(mic.percent_zs0,3);
    [~,~,~, auc(iit,18,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.percent_zs0,3);
    [~,~,~, auc(iit,18,2)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc 
    cc = sum(mean_coh.percent_zs0,3);
    [~,~,~, auc(iit,18,3)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc 
    

end

%% mean and std of AUC 
for ip = 1:18
    for im = 1:3
        median_auc(ip,im) = median(auc(:,ip,im));
        std_auc(ip,im) = std(auc(:,ip,im)); 
        [p_auc(ip,im),h_auc(ip,im),stats] = signrank(auc(:,ip,im),0.5,'tail','right');
        t_auc(ip,im) = stats.signedrank;
    end
end

figure
bar(median_auc)
xTicks = 1:18;
xticklabels = {'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline','max corr','90 percent corr', '1 zs0',...
    '2 zs0', '3 zs0', '4 zs0', '5 zs0', 'max zso', 'percent zs0'};
ylabel('auc')
set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
figure
bar(-log10(p_auc))
set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
ylabel('-log10(p)')

%% corrs  

corrs(:,8:9,3:4)=nan;
for ip = 1:11
    for im = 1:4
        median_corrs(ip,im) = nanmedian(corrs(:,ip,im));
        
    end
end

bar(median_corrs)
xTicks = 1:9;
xticklabels = {'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline'};
ylabel('Pearson correlations')
set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
legend('voxmim','voxmic','voxmeancoh','voxnpcs')

%% variance explained 

for ipip = 1:5 
    median_varex(ipip) = median(varex(:,ipip));
end 

plot(median_varex)
%% save
clear cc
outname= sprintf('/home/bbci/data/haufe/Franziska/data/mimsim_ana_snr0%d_30nit.mat',round(isnr*10));
save(outname,'-v7.3')

%% plots

% boxplots
figure
subplot(2,1,1)
boxplot(squeeze(perf_mim_corr'),'Labels',{'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline','max corr','90 percent corr', '1 zs0',...
    '2 zs0', '3 zs0', '4 zs0', '5 zs0', 'max zso', 'percent zs0'})
grid on
title('MIM correlation with ground truth')
xlabel('Pipelines')
ylabel('correlation coefficient')

% subplot(2,1,2)
% boxplot(squeeze(perf_mim_corr_max'),'Labels',{'1 pc','2 pc', '3 pc',...
%     '4 pc','5 pc','max','90 percent','case2','baseline'})
% grid on
% title('maxima of MIM correlation with ground truth')
% xlabel('Pipelines')
% ylabel('correlation coefficient')

% figure
subplot(2,1,2)
boxplot(squeeze(perf_mic_corr'),'Labels',{'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline'})
grid on
title('MIC correlation with ground truth')
xlabel('Pipelines')
ylabel('correlation coefficient')

% subplot(2,1,2)
% boxplot(squeeze(perf_mic_corr_max'),'Labels',{'1 pc','2 pc', '3 pc',...
%     '4 pc','5 pc','max','90 percent','case2','baseline'})
% grid on
% title('maxima of MIC correlation with ground truth')
% xlabel('Pipelines')
% ylabel('correlation coefficient')

%% hist
figure; 
subplot(2,1,1)
hist(perf_mim_corr_max')
title('maxima of MIC correlation with gt for 9 pipelines')
xlabel('correlation coefficient')
subplot(2,1,2)
hist(perf_mic_corr_max')
title('maxima of MIM correlation with gt for 9 pipelines')
xlabel('correlation coefficient')

