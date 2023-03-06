function fp_figure_S2
% Plot supplementary figure 2

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

DIRDATA = './mim_sim5/';
DIRFIG = './figures/mimsim_ana/mim_sim5/Manuscript/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

labs = {'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};
%%
clearvars -except iname name DIRDATA DIRFIG labs

%default paramenters
nit = 100;
iInt = 2;
iReg=1;
isnr=0.6;
iss = 0.5;
ilag=2;
ifilt='l';
dimred = 'p';
np = 6;

%% load
for iit= 1:nit
    
    inname = sprintf('pr_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_iter%d'...
        ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,dimred, iit);
    
    load([DIRDATA inname '.mat'])
    
    PR{1}(iit,:) = pr_mim;
    PR{2}(iit,:) = pr_mic;
    PR{3}(iit,:) = pr_aCoh;
    PR{4}(iit,:) = pr_iCoh;
    PR{5}(iit,:) = pr_absgc;
    PR{6}(iit,:) = pr_posgc;
    PR{7}(iit,:) = pr_abstrgc;
    PR{8}(iit,:) = pr_postrgc;
    
end

%% plot

figure
figone(50,28)
ipip1=1;
ipip2=1;

pips = [1:8 21 22 10 9]; %fixed, sumvox, summim, central, baseline
npips = 12;

for icon = 1:length(PR) %MIM, MIC, aCoh, iCoh, absgc,posgc,abstrgc, postrgc
    
    for ipip = pips
        
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
        
        if icon >2 && ipip == 10
            
            xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                '99%','Mean+FC','Central','FC+mean','TrueVox'};
            
            if icon == length(PR)
                subplot(length(PR),npips,ipip1)
                view([-90 -90]);
                set(gca, 'Xdir', 'reverse');
                set(gca, 'XLim', [0 1]);
                ylim([-0.75 2])
                htit = title(xtitles(ipip2));
                htit.Position(1) = -0.21;
                set(gca,'ytick',[])
                box off
                set(gca,'xticklabel',{[]})
                set(gca,'XColor','none','YColor','none','TickDir','out')
                set(gca,'Clipping','Off')
                for ix = xt
                    hu = line([ix ix],[2.2 -0.75]);
                    set(hu, 'color',[0.9 0.9 0.9])
                    uistack(hu,'bottom')
                end
                hu = line([0 0],[2.2 -0.75]);
                set(hu, 'color',[0 0 0])
                ipip2=ipip2+1;
            end
            
            ipip1=ipip1+1;
        else
            data1 = PR{icon}(:,ipip);
            mean_pr(icon,ipip) = mean(data1);
            imlab = 'PR';
            imlab1 = 'PR';
            
            %select colors
            if ipip<=np
                cl = [0.8 0.7 0.6];
            elseif ipip <np+3 || ipip == 11 || ipip == 12
                cl = [0.8 0.4 0.5];
            elseif ipip == 21
                cl = [0.4 0.6 0.7];
            elseif ipip == 22
                cl = [0.5 0.7 0.5];
            elseif ipip==9
                cl = [0.8 0.8 0.8];
            elseif ipip == 10
                cl =[0.4 0.5 0.6];
            end
            
            subplot(length(PR),npips,ipip1)
            
            [h, u] = fp_raincloud_plot_a(data1, cl, 1,0.2, 'ks');
            view([-90 -90]);
            set(gca, 'Xdir', 'reverse');
            set(gca, 'XLim', [0 1]);
            ylim([-0.75 2])
            
            xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                '99%','Mean+FC','Central','FC+mean','TrueVox'};
            
            if icon == length(PR)
                htit = title(xtitles(ipip2));
                htit.Position(1) = -0.21;
                ipip2=ipip2+1;
            end
            set(gca,'ytick',[])
            box off
            
            if ipip==1
                xlabel([labs{icon} ' ' imlab1])
                
                set(gca,'Clipping','Off')
                xt = xticks;
                for ix = xt
                    hu = line([ix ix],[2 -40]);
                    
                    set(hu, 'color',[0.9 0.9 0.9])
                    uistack(hu,'bottom')
                end
                hu1 = line([0 0],[2 -40]);
                set(hu1, 'color',[0 0 0])
                
            else
                set(gca,'xticklabel',{[]})
                set(gca,'XColor','none','YColor','none','TickDir','out')
                set(gca,'Clipping','Off')
                for ix = xt
                    hu = line([ix ix],[2.2 -0.75]);
                    set(hu, 'color',[0.9 0.9 0.9])
                    uistack(hu,'bottom')
                end
                hu = line([0 0],[2.2 -0.75]);
                set(hu, 'color',[0 0 0])
            end
            
            ipip1=ipip1+1;
            
        end
    end
end












