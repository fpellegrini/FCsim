function fp_figure3
% Plot different FC metrics 

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

DIRDATA = './mim_sim5/';
DIRFIG = './figures/mimsim_ana/mim_sim5/Manuscript/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

labs = {'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};

%default paramenters
nit = 100;
iInt = 2;
iReg=1;
isnr=0.6;
iss = 0.5;
ilag=2;
ifilt='l';
dimred = 'p';
ipip = 3; %only for fixPC pipeline with 3 PCs

%% load data
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
figone(8,24)
oo = 1;
for icon = [3 4 1 2 5 6 7 8 ] %Coherence, iCOH, MIM, MIC, absgc,posgc, abstrgc, postrgc
    
    data1 = PR{icon}(:,ipip);
    mean_pr(oo) = mean(data1);
    pr_all(oo,:) = data1; 
    imlab1 = 'PR';
    
    cl = [0.8 0.7 0.6];
        
    subplot(1,8,oo)
    
    [h, u] = fp_raincloud_plot_a(data1, cl, 1,0.2, 'ks');
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    set(gca, 'XLim', [0 1]);
    
    htit = title(labs{icon});
    htit.Position(1) = -0.12;
    set(gca,'ytick',[])
    ylim([-0.75 2])
    xlim([0 1])
    box off
    
    if oo==1
        xlabel([imlab1])
        set(gca,'Clipping','Off')
        xt = xticks;
        for ix = xt
            hu = line([ix ix],[2 -25]);
            set(hu, 'color',[0.9 0.9 0.9])
            uistack(hu,'bottom')
        end
        hu1 = line([0 0],[2 -25]);
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
    oo=oo+1;
        
end
%% save

outname = [DIRFIG 'figure3.eps'];
print(outname,'-depsc');
close all

%% test performance differences

%Coherence, iCOH, MIM, MIC, absgc,posgc, abstrgc, postrgc

%coherence vs mim 
d_coh = pr_all(1,:); 
d_mim = pr_all(3,:); 
[p1,~,stats] = signrank(d_coh,d_mim,'tail','left')
t1 = sign(stats.zval);

d_gcdet = pr_all(5,:); 
[p2,~,stats] = signrank(d_gcdet,d_mim,'tail','left')
t2 = sign(stats.zval);


d_gcdir = pr_all(6,:);
d_trgcdir = pr_all(8,:);
[p3,~,stats] = signrank(d_gcdir,d_trgcdir,'tail','left')
t3 = sign(stats.zval);








