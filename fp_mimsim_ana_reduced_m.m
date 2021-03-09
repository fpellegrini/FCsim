function fp_mimsim_ana_reduced_m

DIRIN = '~/data/mim_sim3/';


name = {...
    'ip1_default';...
    'ip2_int2';...
    'ip2_int3';...
    'ip2_int4';...
    'ip2_int5';...
    'ip3_iReg2';...
    'ip4_isnr01';...
    'ip4_isnr03';...
    'ip4_isnr07';...
    'ip4_isnr09';...
    'ip5_iss0';...
    'ip5_iss025';...
    'ip5_iss075';...
    'ip5_iss1';...
    'ip6_lag1';...
    'ip7_dics';...
    'ip7_eloreta_reg';...
    'ip8_hemisym1'};


%%

for iname = 1
    
    clearvars -except iname name DIRIN
    
    %default paramenters
    nit = 100;
    iInt = 2;
    iReg=1;
    isnr=0.7;
    iss = 0.5;
    ilag=2;
    ihemi=0;
    ifilt='l';
    
    if iname>2 && iname<6
        iInt = iname;
    else
        switch iname
            case 1
                iInt = 2;
            case 2
                iInt = 1;
            case 6
                iReg = 2;
            case 7
                isnr = 0.1;
            case 8
                isnr = 0.3;
            case 9
                isnr = 0.5;
            case 10
                isnr = 0.9;
            case 11
                iss=0;
            case 12
                iss = 0.25;
            case 13
                iss = 0.75;
            case 14
                iss = 1;
            case 15
                ilag = 1;
            case 16
                ifilt = 'd';
            case 17
                ifilt = 'e';
            case 18
                ihemi = 1;
        end
    end
    
    band = filt.iband;
    np = 6;
    
    its = 1:100;
    its([31,36,39,43,6,62,65,76])=[];
    
    %%
    for iit= its
        inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d_iter%d'...
            ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,ihemi,iit);
        
        load([DIRIN inname '.mat'])
        
        
<<<<<<< HEAD
        for ipip = 1:np
            cc = sum(mic.fixed{ipip}(:,:,band),3);
=======
        for ipip = 1:6
            cc = sum(mic.fixed{ipip},3);
>>>>>>> 3374706c86a57d7bbb72cf2df6112dff552190ab
            [mrr(iit,ipip,1), pr(iit,ipip,1),hk(iit,ipip,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mic.fixed{ipip}(:,:,band),3);
            [mrr(iit,ipip,2),pr(iit,ipip,2), hk(iit,ipip,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_icoh.fixed{ipip}(:,:,band),3);
            [mrr(iit,ipip,3),pr(iit,ipip,3), hk(iit,ipip,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_acoh.fixed{ipip}(:,:,band),3);
            [mrr(iit,ipip,4),pr(iit,ipip,4), hk(iit,ipip,4)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear c
            

            corrs(iit,ipip,1) = to_save.fixed{ipip}.corr_voxmic;
            corrs(iit,ipip,2) = to_save.fixed{ipip}.corr_voxmim;
            corrs(iit,ipip,3) = to_save.fixed{ipip}.corr_voxmeancoh;
            corrs(iit,ipip,4) = to_save.fixed{ipip}.corr_voxmeanabscoh;
            corrs(iit,ipip,5) = corr(to_save.fixed{ipip}.npcs,to_save.nvoxroi');

            varex(iit,ipip) = to_save.fixed{ipip}.var_explained;
            
        end
        
        %
        cc = sum(mic.max(:,:,band),3);
        [mrr(iit,np+1,1),pr(iit,np+1,1), hk(iit,np+1,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mim.max(:,:,band),3);
        [mrr(iit,np+1,2),pr(iit,np+1,2), hk(iit,np+1,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mean_icoh.max(:,:,band),3);
        [mrr(iit,np+1,3),pr(iit,np+1,3), hk(iit,np+1,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mean_acoh.max(:,:,band),3);
        [mrr(iit,np+1,4),pr(iit,np+1,4), hk(iit,np+1,4)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        
        corrs(iit,np+1,1) = to_save.max.corr_voxmic;
        corrs(iit,np+1,2) = to_save.max.corr_voxmim;
        corrs(iit,np+1,3) = to_save.max.corr_voxmeancoh;
        corrs(iit,np+1,4) = to_save.max.corr_voxmeanabscoh;
        corrs(iit,np+1,5) = corr(to_save.max.npcs',to_save.nvoxroi');
        %
        cc = sum(mic.percent(:,:,band),3);
        [mrr(iit,np+2,1), pr(iit,np+2,1),hk(iit,np+2,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mim.percent(:,:,band),3);
        [mrr(iit,np+2,2), pr(iit,np+2,2),hk(iit,np+2,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mean_icoh.percent(:,:,band),3);
        [mrr(iit,np+2,3), pr(iit,np+2,3),hk(iit,np+2,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mean_acoh.percent(:,:,band),3);
        [mrr(iit,np+2,4), pr(iit,np+2,4),hk(iit,np+2,4)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        
        
        corrs(iit,np+2,1) = to_save.percent.corr_voxmic;
        corrs(iit,np+2,2) = to_save.percent.corr_voxmim;
        corrs(iit,np+2,3) = to_save.percent.corr_voxmeancoh;
        corrs(iit,np+2,4) = to_save.percent.corr_voxmeanabscoh;
        corrs(iit,np+2,5) = corr(to_save.percent.npcs',to_save.nvoxroi');
        
        cc = sum(mic.baseline(:,:,band),3);
        [mrr(iit,np+4,1), pr(iit,np+4,1),hk(iit,np+4,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        cc = sum(mim.baseline(:,:,band),3);
        [mrr(iit,np+4,2),pr(iit,np+4,2), hk(iit,np+4,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
        clear cc
        
        corrs(iit,np+4,1) = to_save.baseline.corr_voxmic;
        corrs(iit,np+4,2) = to_save.baseline.corr_voxmim;
        
        if iname == 1
            cc= sum(mic.case2(:,:,band),3);
            [mrr(iit,np+3,1),pr(iit,np+3,1), hk(iit,np+3,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc= sum(mim.case2(:,:,band),3);
            [mrr(iit,np+3,2),pr(iit,np+3,2), hk(iit,np+3,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            
            %corrected
            cc = sum(mic.max_corrected(:,:,band),3);
            [mrr(iit,np+5,1),pr(iit,np+5,1), hk(iit,np+5,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mim.max_corrected(:,:,band),3);
            [mrr(iit,np+5,2),pr(iit,np+5,2), hk(iit,np+5,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_icoh.max_corrected(:,:,band),3);
            [mrr(iit,np+5,3),pr(iit,np+5,3), hk(iit,np+5,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_acoh.max_corrected(:,:,band),3);
            [mrr(iit,np+5,4),pr(iit,np+5,4), hk(iit,np+5,4)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            
            corrs(iit,np+5,1) = to_save.max_corrected.corr_voxmic;
            corrs(iit,np+5,2) = to_save.max_corrected.corr_voxmim;
            corrs(iit,np+5,3) = to_save.max_corrected.corr_voxmeancoh;
            corrs(iit,np+5,4) = to_save.max_corrected.corr_voxmeanabscoh;
            %
            cc = sum(mic.percent_corrected(:,:,band),3);
            [mrr(ii,np+6,1), pr(iit,np+6,1),hk(iit,np+6,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mim.percent_corrected(:,:,band),3);
            [mrr(iit,np+6,2), pr(iit,np+6,2),hk(iit,np+6,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_icoh.percent_corrected(:,:,band),3);
            [mrr(iit,np+6,3), pr(iit,np+6,3),hk(iit,np+6,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_acoh.percent_corrected(:,:,band),3);
            [mrr(iit,np+6,4), pr(iit,np+6,4),hk(iit,np+6,4)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
        
        
            corrs(iit,np+6,1) = to_save.percent_corrected.corr_voxmic;
            corrs(iit,np+6,2) = to_save.percent_corrected.corr_voxmim;
            corrs(iit,np+6,3) = to_save.percent_corrected.corr_voxmeancoh;
            corrs(iit,np+6,4) = to_save.percent_corrected.corr_voxmeanabscoh;
            
            
            %zs0
            
<<<<<<< HEAD
            for ipip = 1:np
                cc = sum(mic.fixed_zs0{ipip}(:,:,band),3);
                [mrr(iit,np+6+ipip,1), pr(iit,np+6+ipip,1),hk(iit,np+6+ipip,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
=======
            for ipip = 1:6
                cc = sum(mic.fixed_zs0{ipip},3);
                [mrr(iit,11+ipip,1), pr(iit,11+ipip,1),hk(iit,11+ipip,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
>>>>>>> 3374706c86a57d7bbb72cf2df6112dff552190ab
                clear cc
                cc = sum(mic.fixed_zs0{ipip}(:,:,band),3);
                [mrr(iit,np+6+ipip,2),pr(iit,np+6+ipip,2), hk(iit,np+6+ipip,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
                clear cc
                cc = sum(mean_icoh.fixed_zs0{ipip}(:,:,band),3);
                [mrr(iit,np+6+ipip,3),pr(iit,np+6+ipip,3), hk(iit,np+6+ipip,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
                clear cc
                cc = sum(mean_acoh.fixed_zs0{ipip}(:,:,band),3);
                [mrr(iit,np+6+ipip,4),pr(iit,np+6+ipip,4), hk(iit,np+6+ipip,4)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
                clear c
                
            end

            %
            cc = sum(mic.max_zs0(:,:,band),3);
            [mrr(iit,2*np+7,1),pr(iit,2*np+7,1), hk(iit,2*np+7,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mim.max_zs0(:,:,band),3);
            [mrr(iit,2*np+7,2),pr(iit,2*np+7,2), hk(iit,2*np+7,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_icoh.max_zs0(:,:,band),3);
            [mrr(iit,2*np+7,3),pr(iit,2*np+7,3), hk(iit,2*np+7,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_acoh.max_zs0(:,:,band),3);
            [mrr(iit,2*np+7,4),pr(iit,2*np+7,4), hk(iit,2*np+7,4)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            
            %
            cc = sum(mic.percent_zs0(:,:,band),3);
            [mrr(iit,2*np+8,1), pr(iit,2*np+8,1),hk(iit,2*np+8,1)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mim.percent_zs0(:,:,band),3);
            [mrr(iit,2*np+8,2), pr(iit,2*np+8,2),hk(iit,2*np+8,2)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_icoh.percent_zs0(:,:,band),3);
            [mrr(iit,2*np+8,3), pr(iit,2*np+8,3),hk(iit,2*np+8,3)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            cc = sum(mean_acoh.percent_zs0(:,:,band),3);
            [mrr(iit,2*np+8,4), pr(iit,2*np+8,4),hk(iit,2*np+8,4)] = fp_mrr_hk(cc, iroi_seed,iroi_tar);
            clear cc
            
        end
        
    end
    
    %%
    if iname > 1
        mrr(:,np+3,:)=[];
        pr(:,np+3,:)=[];
        hk(:,np+3,:)=[];
    end
    
    %%
    mrr([31,36,39,43,6,62,65,76],:,:)=[];
    pr([31,36,39,43,6,62,65,76],:,:)=[];
    hk([31,36,39,43,6,62,65,76],:,:)=[];
    %%
    outname = [DIRIN 'mrr_mim3_' name{iname} '.mat'];
    save(outname,'mrr','pr','hk','corrs','varex','-v7.3')
end
%%
%
% a = {'1zs','2zs', '3zs','4zs','5zs',...
%     '99% zs','90% zs','sumVox','baseline'};
%
% b = {'mic','mim'};
%
%
% for jj = 1:2
%     figone(30,50)
%     for ii = 1:9
%         subplot(3,3,ii)
%         hist(pr(:,ii,jj))
%         title(['mrr ' a{ii} ' ' b{jj}])
%         ylim([0 100])
%
%     end
% %     saveas(gcf,['../figures/mimsim_ana/mrrs_' b{jj}],'png')
% %     close all
% end
%
% for jj = 1:2
%     figone(30,50)
%     for ii = 1:9
%         subplot(3,3,ii)
%         hist(hk(:,ii,jj))
%         title(['hits at 2 ' a{ii} ' ' b{jj}])
%         ylim([0 100])
%         
%     end
% %     saveas(gcf,['../figures/mimsim_ana/hk_' b{jj}],'png')
% %     close all
% end
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
