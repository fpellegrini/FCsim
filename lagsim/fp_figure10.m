DIRDATA = './mim_sim5_lag/';
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
    'ip7_dics'};

labs = {'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};

im = 2; %measures: MRR, PR, EM3
ipip = 3; %only third pipeline was calculated
%%
icon = 8; %labs = {'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};
o=1;
figure
figone(8,24)
iname = 1;

clearvars -except iname name DIRDATA DIRFIG labs o im icon ipip mean_pr

%default paramenters
nit = 100;
iInt = 2;
iReg=1;
isnr=0.6;
iss = 0.5;
ifilt='l';
dimred='p';

np = 6;
a=[];
%%
for ilag = 1:6
    for iit= 1:nit
        
        %% load data
        if ilag == 6 
            inname = sprintf('pr_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_iter%d'...
                ,iInt,iReg,isnr*10,iss*10, 2,ifilt,dimred, iit);
            load(['./mim_sim5/' inname '.mat'])
        else
            inname = sprintf('pr_iInt%d_iReg%d_snr0%d_iss0%d_filt%s_iter%d_lag%d'...
                ,iInt,iReg,isnr*10,iss*10, ifilt, iit, ilag);
            load([DIRDATA inname '.mat'])
        end               
        
        PR{1}(iit,:) = pr_mim(1:3);
        PR{2}(iit,:) = pr_mic(1:3);
        PR{3}(iit,:) = pr_aCoh(1:3);
        PR{4}(iit,:) = pr_iCoh(1:3);
        PR{5}(iit,:) = pr_absgc(1:3);
        PR{6}(iit,:) = pr_posgc(1:3);
        PR{7}(iit,:) = pr_abstrgc(1:3);
        PR{8}(iit,:) = pr_postrgc(1:3);        
    end

    
    data1 = PR{icon}(:,ipip);
    mean_pr(o) = mean(data1);
    imlab = 'PR';
    imlab1 = 'PR';
    
    %select colors
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
    
    %% plot
    subplot(1,6,o)
    
    [h, u] = fp_raincloud_plot_a(data1, cl, 1,0.2, 'ks');
    
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    set(gca, 'XLim', [0 1]);
    
    titles = {'2 ms' ,'4 ms' ,'6 ms','8 ms','10 ms','50-200 ms'};
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
    
    o=o+1;
    
end

%% save
outname = [DIRFIG imlab '_' labs{icon} '_lags.eps'];
print(outname,'-depsc');
close all



