function fp_figure9

DIRDATA = './mim_sim5/';
DIRFIG = './figures/mimsim_ana/mim_sim5/Manuscript/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

name = {...
    'ip1_default';...
    'ip2_int1';...
    'ip2_int3';...
    'ip2_int4';...
    'ip2_int5';...
    'ip3_iReg2';...
    'ip4_isnr05';...
    'ip4_isnr09';...
    'ip5_iss0';...
    'ip5_iss025';...
    'ip5_iss075';...
    'ip5_iss1';...
    'ip6_lag1';...
    'ip7_eloreta_reg';...%14
    'ip7_champ';...
    'ip7_champ_reg';...
    'ip8_ssd';...
    'ip7_dics';...
    'ip7_che';
    'ip7_cho'};%20

labs = {'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};
[cb] = cbrewer2('spectral', 11);
cb1 = cbrewer2('Set1',9);

iname = 6;
im = 2;
%%
clearvars -except iname name DIRDATA DIRFIG labs cb cb1 mean_pr

%default paramenters
nit = 100;
iInt = 2;
isnr=0.6;
iss = 0.5;
ilag=2;
ihemi=0;
ifilt='l';
iReg = 2;
dimred = 'p';

np = 6;

its = [1:100];
a=[];

%%
for iit= its
    
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

%%


for icon = [1 8] %labs = {'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};
    
    figure
    figone(8,24)
    
    pips =  [1:6 ]; %9
    npips = 6; %7
    pip1=1;
    
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
        
        data1 = PR{icon}(:,ipip);
        mean_pr(pip1,icon) = mean(data1);
        imlab = 'PR';
        imlab1 = 'PR';
        
        
        if ipip<=np
            cl = [0.8 0.7 0.6];
        elseif ipip <np+3 || ipip == 11 || ipip == 12
            cl = [0.8 0.4 0.5];
        elseif ipip == 10
            cl = [0.4 0.6 0.7];
        elseif ipip == 9
            cl = [0.8 0.8 0.8];
        end
        
        subplot(1,npips,pip1)
        
        [h, u] = fp_raincloud_plot_a(data1, cl, 1,0.2, 'ks');
        view([-90 -90]);
        set(gca, 'Xdir', 'reverse');
        set(gca, 'XLim', [0 1]);
        
        xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%', '99%','TrueVox'};
        htit = title(xtitles(ipip));
        htit.Position(1) = -0.12;
        set(gca,'ytick',[])
        ylim([-0.75 2])
        xlim([0 1])
        box off
        
        if pip1==1
            xlabel([imlab1])
            set(gca,'Clipping','Off')
            xt = xticks;
            for ix = xt
                hu = line([ix ix],[2 -18]);
                set(hu, 'color',[0.9 0.9 0.9])
                uistack(hu,'bottom')
            end
            hu1 = line([0 0],[2 -18]);
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
        
        pip1 = pip1+1;
        
    end
    
    ylim([-0.75 2])
    %%
    outname = [DIRFIG name{iname} '_' imlab '_' labs{icon} '.eps'];
    print(outname,'-depsc');
    close all
    
    
    
end








