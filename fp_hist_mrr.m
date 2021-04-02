DIRDATA = './mim_sim3/';
DIRFIG = './figures/mimsim_ana/mim_sim3/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

name = {...
    'ip1_default';...
    'ip2_int1';...
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
    'ip7_champ';...
    'ip7_eloreta_reg'};

labs = {'MIC','MIM','Mean abscoh','mean icoh','absGC','posGC'};
[cb] = cbrewer2('spectral', 11);
cb1 = cbrewer2('Set1',9);

%%

for iname = 1%:numel(name)
    
    
    clearvars -except iname name DIRDATA DIRFIG labs cb cb1
    
    %default paramenters
    nit = 100;
    iInt = 2;
    iReg=1;
    isnr=0.7;
    iss = 0.5;
    ilag=2;
    ihemi=0;
    ifilt='l';
    
    if iname>1 && iname<6
        iInt = iname;
    else
        switch iname
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
                ifilt = 'c';
            case 17
                ifilt = 'e';
        end
    end
    
    
    np = 6;
    
    its = [1:100];
    a=[];
    
    %%
    for iit= its
        
        try
            inname = sprintf('mrr_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
                ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
            
            load([DIRDATA inname '.mat'])
            
            MRR{1}(iit,:) = mrr_mim;
            MRR{2}(iit,:) = mrr_mic;
            MRR{3}(iit,:) = mrr_aCoh;
            MRR{4}(iit,:) = mrr_iCoh;
            MRR{5}(iit,:) = mrr_absgc;
            MRR{6}(iit,:) = mrr_posgc;
            
            PR{1}(iit,:) = pr_mim;
            PR{2}(iit,:) = pr_mic;
            PR{3}(iit,:) = pr_aCoh;
            PR{4}(iit,:) = pr_iCoh;
            PR{5}(iit,:) = pr_absgc;
            PR{6}(iit,:) = pr_posgc;
            
        catch
            a = [a iit];
        end
    end
    
    for ii = 1:6
        MRR{ii}(a,:) = [];
        PR{ii}(a,:) = [];
    end
    
    %%
    
    
    for icon = 1:length(MRR)
        
        
        %%
        figure
        figone(15,50)
        
        if icon < 5
            npips = 12;
        else 
            npips = 9;
        end
        
        for ipip = 1:npips
            
            %1 to 6: fixed
            %7: 90%
            %8: 99%
            %9: baseline
            %10: sumVox
            %11: 90% corrected
            %12: 99% corrected
            %13 to 20: equal to 1 to 8 but with zscoring
            %21: sum, then MIM
            %22: central voxel
            
            
            subplot(1,npips,ipip)
            
            if ipip<=np
                cl = [0.5 0.5 1];%cb(np+1-ipip,:);
            elseif ipip <np+3 || ipip == 11 || ipip == 12
                cl = [1 0.5 0.5];%cb(ipip+2,:);
            elseif ipip == 10
                cl = [0.5 1 0.5];
            elseif ipip == 9
                cl = [0.2 0.2 0.2];
            end
            
            data1 = MRR{icon}(:,ipip);
            [h, u] = fp_raincloud_plot(data1, cl, 1,0.2, 'ks');
            view([-90 -90]);
            set(gca, 'Xdir', 'reverse');
            set(gca, 'XLim', [0 1]);
            
            if iname == 1
                %                     if icon <3
                xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                    '99%','baseline','sumVox','90% corrected','99% corrected'};
                %                     else
                %                         xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                %                             '99%'};
                %                     end
            else
                xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%', '99%','baseline'};
            end
            
            title(xtitles(ipip))
            if ipip~=1
                set(gca,'xtick',[])
            end
            
            xlabel([labs{icon} ' MRR'])
            
        end
        
        %                 xlim([0 103])
        ylim([-0.75 2])
        
        outname = [DIRFIG name{iname} '_' labs{icon}];
        saveas(gcf,outname, 'png')
        close all
        
        
    end
    
    
    if iname ==1
        
        for icon = 1:length(MRR)
            
            
            %%
            figure
            figone(15,50)
            
            if icon <5
                npips = size(MRR{1},2)-12;
            else 
                npips = size(MRR{1},2)-14;
            end
            
            for ipip = 1:npips
                
                
                subplot(1,npips,ipip)
                
                if ipip<=np
                    cl = [0.5 0.5 1];%cb(np+1-ipip,:);
                elseif ipip <np+3 || ipip == 11 || ipip == 12
                    cl = [1 0.5 0.5];%cb(ipip+2,:);
                elseif ipip == 10
                    cl = [0.8 0.8 0.5];
                elseif ipip == 9
                    cl = [0.6 0.2 0.8];
                end
                
                data1 = MRR{icon}(:,ipip+12);
                [h, u] = fp_raincloud_plot(data1, cl, 1,0.2, 'ks');
                view([-90 -90]);
                set(gca, 'Xdir', 'reverse');
                set(gca, 'XLim', [0 1]);
                
                
                %                     if icon <3
                xtitles = {'1PC-zs','2PCs-zs', '3PCs-zs','4PCs-zs','5PCs-zs','6PCs-zs','90%-zs',...
                    '99%-zs','sum+MIM','central voxel'};
                %                     else
                %                         xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                %                             '99%'};
                %                     end
                
                title(xtitles(ipip))
                if ipip~=1
                    set(gca,'xtick',[])
                end
                
                
                xlabel([labs{icon} ' MRR'])
                
            end
            
            %                 xlim([0 103])
            ylim([-0.75 2])
            
            outname = [DIRFIG name{iname} '_2_' labs{icon}];
            saveas(gcf,outname, 'png')
            close all
            
            
        end
        
    end
    
   
    
    
end




    
    

