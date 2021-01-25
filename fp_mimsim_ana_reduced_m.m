function fp_mimsim_ana_reduced_m

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim2/';

%default paramenters
nit = 100;
iInt = 1; 
iReg=1; 
isnr=0.5;
iss = 0.5;
ilag=2;
ihemi=0;
ifilt='d';


%%
for iit= [1:nit]
    inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d_iter%d'...
        ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,ihemi,iit);
    
    load([DIRIN inname '.mat'])
   
    
    for ipip = 1:5
        cc = sum(mic.fixed{ipip},3);
        [mrr(iit,ipip,1), mrrs(iit,ipip,1),hk(iit,ipip,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mic.fixed{ipip},3);
        [mrr(iit,ipip,2),mrrs(iit,ipip,2), hk(iit,ipip,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mean_icoh.fixed{ipip},3);
        [mrr(iit,ipip,3),mrrs(iit,ipip,3), hk(iit,ipip,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
%         
%         corrs(iit,ipip,1) = to_save.fixed{ipip}.corr_voxmic;
%         corrs(iit,ipip,2) = to_save.fixed{ipip}.corr_voxmim; 
%         corrs(iit,ipip,3) = to_save.fixed{ipip}.corr_voxmeancoh;
%         corrs(iit,ipip,4) = to_save.fixed{ipip}.corr_voxmeanabscoh;
%         corrs(iit,ipip,5) = corr(to_save.fixed{ipip}.npcs,to_save.nvoxroi');
%         
%         varex(iit,ipip) = to_save.fixed{ipip}.var_explained;
        
    end

%
    cc = sum(mic.max,3);
    [mrr(iit,6,1),mrrs(iit,6,1), hk(iit,6,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
    clear cc
    cc = sum(mim.max,3);
    [mrr(iit,6,2),mrrs(iit,6,2), hk(iit,6,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
    clear cc
    cc = sum(mean_icoh.max,3);
    [mrr(iit,6,3),mrrs(iit,6,3), hk(iit,6,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
    clear cc
    
%     corrs(iit,6,1) = to_save.max.corr_voxmic;
%     corrs(iit,6,2) = to_save.max.corr_voxmim; 
%     corrs(iit,6,3) = to_save.max.corr_voxmeancoh;
%     corrs(iit,6,4) = to_save.max.corr_voxmeanabscoh;
%     corrs(iit,6,5) = corr(to_save.max.npcs',to_save.nvoxroi');
%     
    cc = sum(mic.percent,3);
    [mrr(iit,7,1), mrrs(iit,7,1),hk(iit,7,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
    clear cc
    cc = sum(mim.percent,3);
    [mrr(iit,7,2), mrrs(iit,7,2),hk(iit,7,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
    clear cc 
    cc = sum(mean_icoh.percent,3);
    [mrr(iit,7,3), mrrs(iit,7,3),hk(iit,7,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
    clear cc 
%     corrs(iit,7,1) = to_save.percent.corr_voxmic;
%     corrs(iit,7,2) = to_save.percent.corr_voxmim; 
%     corrs(iit,7,3) = to_save.percent.corr_voxmeancoh;
%     corrs(iit,7,4) = to_save.percent.corr_voxmeanabscoh;
%     corrs(iit,7,5) = corr(to_save.percent.npcs',to_save.nvoxroi');
  
    cc = sum(mic.baseline,3);
    [mrr(iit,9,1), mrrs(iit,9,1),hk(iit,9,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
    clear cc
    cc = sum(mim.baseline,3);
    [mrr(iit,9,2),mrrs(iit,9,2), hk(iit,9,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
    clear cc
    
%     corrs(iit,8,1) = to_save.baseline.corr_voxmic;
%     corrs(iit,8,2) = to_save.baseline.corr_voxmim; 

%     cc= sum(mic.case2,3);
%      [mrr(iit,8,1),mrrs(iit,8,1), hk(iit,8,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
%      clear cc
%      cc= sum(mim.case2,3);
%      [mrr(iit,8,2),mrrs(iit,8,2), hk(iit,8,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
%      clear cc

end
    
%%
mrr(:,8,:)=[];
mrrs(:,8,:)=[];
hk(:,8,:)=[];
%%
outname = [DIRIN 'mrr_mim2_dics.mat'];
save(outname,'mrr','mrrs','hk','-v7.3')
%%

a = {'1zs','2zs', '3zs','4zs','5zs',...
    '99% zs','90% zs','sumVox','baseline'};

b = {'mic','mim'};


for jj = 1:2
    figone(30,50)
    for ii = 1:9
        subplot(3,3,ii)
        hist(mrrs(:,ii,jj))
        title(['mrr ' a{ii} ' ' b{jj}])
        ylim([0 100])
        
    end
%     saveas(gcf,['../figures/mimsim_ana/mrrs_' b{jj}],'png')
%     close all
end

for jj = 1:2
    figone(30,50)
    for ii = 1:9
        subplot(3,3,ii)
        hist(hk(:,ii,jj))
        title(['hits at 2 ' a{ii} ' ' b{jj}])
        ylim([0 100])
        
    end
%     saveas(gcf,['../figures/mimsim_ana/hk_' b{jj}],'png')
%     close all
end
%     
% %     cc = sum(mic.max_corrected,3);
% %     [~,~,~, auc(iit,10,1)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc
% %     cc = sum(mim.max_corrected,3);
% %     [~,~,~, auc(iit,10,2)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc
% %     cc = sum(mean_icoh.max_corrected,3);
% %     [~,~,~, auc(iit,10,3)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc
% %     corrs(iit,10,1) = to_save.max_corrected.corr_voxmic;
% %     corrs(iit,10,2) = to_save.max_corrected.corr_voxmim; 
% %     corrs(iit,10,3) = to_save.max_corrected.corr_voxmeancoh;
% %     corrs(iit,10,4) = corr(to_save.max.npcs',to_save.nvoxroi');
% %     
% %     cc = sum(mic.percent_corrected,3);
% %     [~,~,~, auc(iit,11,1)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc
% %     cc = sum(mim.percent_corrected,3);
% %     [~,~,~, auc(iit,11,2)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc 
% %     cc = sum(mean_icoh.percent_corrected,3);
% %     [~,~,~, auc(iit,11,3)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc 
% %     corrs(iit,11,1) = to_save.percent_corrected.corr_voxmic;
% %     corrs(iit,11,2) = to_save.percent_corrected.corr_voxmim; 
% %     corrs(iit,11,3) = to_save.percent_corrected.corr_voxmeancoh;
% %     corrs(iit,11,4) = corr(to_save.percent.npcs',to_save.nvoxroi');
%     
%     
% %     for ipip = 1:5
% %         cc = sum(mic.fixed_zs0{ipip},3);
% %         [~, ~,~, auc(iit,11+ipip,1)] = ...
% %             perfcurve(label(:),cc(:),1);
% %         clear cc
% %         cc = sum(mim.fixed_zs0{ipip},3);
% %         [~,~,~, auc(iit,11+ipip,2)] = ...
% %             perfcurve(label(:),cc(:),1);
% %         clear cc
% %         cc = sum(mean_icoh.fixed_zs0{ipip},3);
% %         [~,~,~, auc(iit,11+ipip,3)] = ...
% %             perfcurve(label(:),cc(:),1);
% %         clear cc
% %         
% %        
% %     end
% %      
% %     cc = sum(mic.max_zs0,3);
% %     [~,~,~, auc(iit,17,1)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc
% %     cc = sum(mim.max_zs0,3);
% %     [~,~,~, auc(iit,17,2)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc
% %     cc = sum(mean_icoh.max_zs0,3);
% %     [~,~,~, auc(iit,17,3)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc
% %     
% %     cc = sum(mic.percent_zs0,3);
% %     [~,~,~, auc(iit,18,1)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc
% %     cc = sum(mim.percent_zs0,3);
% %     [~,~,~, auc(iit,18,2)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc 
% %     cc = sum(mean_icoh.percent_zs0,3);
% %     [~,~,~, auc(iit,18,3)] = ...
% %         perfcurve(label(:),cc(:),1);
% %     clear cc 
%     
% 
% end
% 
% %% save auc 
% outname = '/home/bbci/data/haufe/Franziska/data/mim_sim1_auc_ip8_hemisym1';
% save(outname,'auc','-v7.3')
% 
% %% mean and std of AUC 
% for ip = 1:8
%     for im = 1:4
%         median_auc(ip,im) = mean(auc(:,ip,im));
%         std_auc(ip,im) = std(auc(:,ip,im)); 
%         [p_auc(ip,im),h_auc(ip,im),stats] = signrank(auc(:,ip,im),0.5,'tail','right');
%         t_auc(ip,im) = stats.signedrank;
%     end
% end
% 
% figure
% bar(median_auc)
% xTicks = 1:8;
% xticklabels = {'1 pc','2 pc', '3 pc','4 pc','5 pc',...
%     '99%','90%','baseline'};
% ylabel('auc')
% set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
% legend('mic','mim','mean icoh', 'mean abscoh')
% grid on 
% title('same number of PCs in regions of both hemispheres')
% % subplot(2,1,2)
% % bar(-log10(p_auc))
% % set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
% % ylabel('-log10(p)')
% % legend('mic','mim','mean coh')
% % grid on 
% 
% %% corrs  
% corrs(:,8,3:5)=nan;
% %%
% for ip = 1:8
%     for im = 1:5
%         median_corrs(ip,im) = nanmedian(corrs(:,ip,im));
%         
%     end
% end
% 
% bar(median_corrs)
% xTicks = 1:11;
% xticklabels = {'1 pc','2 pc', '3 pc','4 pc','5 pc',...
%     '99%','90%','baseline'};
% 
% ylabel('Pearson correlations')
% set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
% legend('voxmim','voxmic','voxmeanicoh','voxmean abscoh','voxnpcs')
% title('3 interactions')
% 
% %% variance explained
% %%
% for ipip = 1:5 
%     median_varex(ipip) = median(varex(:,ipip));
% end 
% 
% plot(median_varex)
% grid on
% xlabel('Number of fixed PCs')
% ylabel('Variance explained')
% title('variance explained, median across iterations, only brain noise ')
% % %% save
% % clear cc
% % outname= sprintf('/home/bbci/data/haufe/Franziska/data/mimsim_ana1_snr0%d.mat',round(isnr*10));
% % save(outname,'-v7.3')