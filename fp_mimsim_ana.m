function fp_mimsim_ana

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim/';

%default sparamenters
nit = 100;
iInt = 1; 
iReg=1; 
isnr=0.5;
iss = 0.5;
ilag=2;
ihemi=0;
ifilt='l';

%%
for iit=1:nit
    inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d_iter%d'...
        ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,ihemi,iit);
    
    load([DIRIN inname '.mat'])
        
    %ground truth
    gt.mic = abs(imag(gt.mic));
    gt.mim = abs(imag(gt.mim)); 
    
    %there was some bug, so we recalculate PERFORMANCE and BASELINE
    %PERFORMANCE has the dimensions: 2 (mim/mic) x 8 (pipelines) x
    %2 (perfomance measure)
    clear PERFORMANCE BASELINE
    [PERFORMANCE, BASELINE] = fp_get_performance(gt, mic, mim, params);

    perf_mim_corr(:,iit) = [squeeze(PERFORMANCE(1,:,1)) squeeze(BASELINE(1,1))];
    perf_mic_corr(:,iit) = [squeeze(PERFORMANCE(2,:,1)) squeeze(BASELINE(1,1))];
%     perf_mim_corr_max(:,iit) = [squeeze(PERFORMANCE(1,:,2)) squeeze(BASELINE(1,1))];
%     perf_mic_corr_max(:,iit) = [squeeze(PERFORMANCE(2,:,2)) squeeze(BASELINE(1,1))];
    
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
    end
    cc = sum(mic.max,3);
    [~,~,~, auc(iit,6,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.max,3);
    [~,~,~, auc(iit,6,2)] = ...
        perfcurve(label(:),cc(:),1);
    cc = sum(mic.percent,3);
    [~,~,~, auc(iit,7,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.percent,3);
    [~,~,~, auc(iit,7,2)] = ...
        perfcurve(label(:),cc(:),1);
    cc = sum(mic.case2,3);
    [~,~,~, auc(iit,8,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.case2,3);
    [~,~,~, auc(iit,8,2)] = ...
        perfcurve(label(:),cc(:),1);
    cc = sum(mic.baseline,3);
    [~,~,~, auc(iit,9,1)] = ...
        perfcurve(label(:),cc(:),1);
    clear cc
    cc = sum(mim.baseline,3);
    [~,~,~, auc(iit,9,2)] = ...
        perfcurve(label(:),cc(:),1);
          

end

%%
for ip = 1:9
    for im = 1:2 
        median_auc(ip,im) = median(auc(:,ip,im));
        std_auc(ip,im) = std(auc(:,ip,im)); 
        [p_auc(ip,im),h_auc(ip,im),stats] = signrank(auc(:,ip,im),0.5,'tail','right');
        t_auc(ip,im) = stats.signedrank;
    end
end

figure
bar(median_auc)
title('median auc')
xlabel('pipelines')
ylabel('auc')
legend('mic','mim')
figure
bar(-log10(p_auc))
title('-log10(p auc)')
xlabel('pipelines')
ylabel('-log10(p auc)')
legend('mic','mim')
%% save
clear cc
outname= '/home/bbci/data/haufe/Franziska/data/mimsim_ana.mat';
save(outname,'-v7.3')
%%
%% plots

% boxplots
figure
subplot(2,1,1)
boxplot(squeeze(perf_mim_corr'),'Labels',{'1 pc','2 pc', '3 pc',...
    '4 pc','5 pc','max','90 percent','case2','baseline'})
grid on
title('MIM correlation with ground truth')
xlabel('Pipelines')
ylabel('correlation coefficient')
% 
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
% 
% subplot(2,1,2)
% boxplot(squeeze(perf_mic_corr_max'),'Labels',{'1 pc','2 pc', '3 pc',...
%     '4 pc','5 pc','max','90 percent','case2','baseline'})
% grid on
% title('maxima of MIC correlation with ground truth')
% xlabel('Pipelines')
% ylabel('correlation coefficient')

%%
%hist
figure; 
subplot(2,1,1)
hist(perf_mim_corr_max')
title('maxima of MIC correlation with gt for 9 pipelines')
xlabel('correlation coefficient')
subplot(2,1,2)
hist(perf_mic_corr_max')
title('maxima of MIM correlation with gt for 9 pipelines')
xlabel('correlation coefficient')

