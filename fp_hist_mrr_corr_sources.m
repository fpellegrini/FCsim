function fp_hist_mrr_corr_sources

DIRDATA = './mim_sim5/';
DIRFIG = './figures/mimsim_ana/mim_sim5/Manuscript/';
if ~exist(DIRFIG); mkdir(DIRFIG); end


labs = {'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};
titles = {'eLORETA','LCMV','DICS','Champagne'};

ipip = 3;
%%
for icon = 1:8 %{'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};
    o=1;
    figure
    figone(8,18)
    for iname = [14 1]
        
        clearvars -except iname name DIRDATA DIRFIG labs o im icon ipip titles xt mean_pr
        
        %default paramenters
        nit = 100;
        iInt = 2;
        iReg=1;
        isnr=0.6;
        iss = 0.5;
        ilag=2;
        ifilt='l';
        dimred = 'p';
        
        if iname==2
            iInt = 1;
        elseif iname>2 && iname<6
            iInt = iname;
        else
            switch iname
                case 6
                    iReg = 2;
                case 7
                    isnr = 0.5;
                case 8
                    isnr = 0.9;
                case 9
                    iss=0;
                case 10
                    iss = 0.25;
                case 11
                    iss = 0.75;
                case 12
                    iss = 1;
                case 13
                    ilag = 1;
                case 14
                    ifilt = 'e';
                case 15
                    ifilt = 'c';
                case 16
                    ifilt = 'cr';
                case 18
                    ifilt = 'd';
                    dimred = 'p';
                case 17
                    dimred = 's';
            end
        end
        
        
        np = 6;
        
        its = [1:nit];
        a=[];
        
        %%
        for iit= its
            
            
            inname = sprintf('pr_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_corr_iter%d'...
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
        data1 = PR{icon}(:,ipip);
        mean_pr(o) = mean(data1);
        imlab = 'PR';
        imlab1 = 'PR';
        
        
        if ipip<=np
            cl = [0.8 0.7 0.6];
        elseif ipip <np+3 || ipip == 11 || ipip == 12
            cl = [0.8 0.4 0.5];
        elseif ipip == 10 || ipip == 21
            cl = [0.4 0.6 0.7];
        elseif ipip == 22
            cl = [0.5 0.7 0.5];
        elseif ipip==9
            cl = [0.8 0.8 0.8];
        end
        
        subplot(1,4,o)
        
        [h, u] = fp_raincloud_plot_a(data1, cl, 1,0.2, 'ks');
        view([-90 -90]);
        set(gca, 'Xdir', 'reverse');
        set(gca, 'XLim', [0 1]);
        
        htit = title(titles{o});
        htit.Position(1) = -0.12;
        set(gca,'ytick',[])
        ylim([-0.75 2])
        box off
        
        if o==1
            xlabel([labs{icon} ' ' imlab1])
            set(gca,'Clipping','Off')
            xt = xticks;
            for ix = xt
                hu = line([ix ix],[2 -10]);
                set(hu, 'color',[0.9 0.9 0.9])
                uistack(hu,'bottom')
            end
            hu1 = line([0 0],[2 -10]);
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
        
        o=o+1;
        
    end
    
    %%
    outname = [DIRFIG imlab '_' labs{icon} '_correlated_sources_3pcs.eps'];
    print(outname,'-depsc');
    
    close all
    
end






