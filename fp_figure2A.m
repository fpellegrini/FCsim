function fp_figure2A

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
    'ip7_eloreta_reg';...
    'ip7_champ';...
    'ip7_champ_reg';...
    'ip8_ssd';...
    'ip7_dics'};

labs = {'MIM','MIC','Coherence','iCOH','GC-det','GC-dir','TRGC-det','TRGC-dir'};
%%
iname = 1;


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

its = [1:nit];
a=[];
im = 2; %measure: PR

%%
for iit= its
    
    try
%         if iname == 18 | iname == 17
        inname = sprintf('pr_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_iter%d'...
                ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,dimred, iit);
%         else
%             inname = sprintf('mrr_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
%                 ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
%         end
        
        load([DIRDATA inname '.mat'])
        
        PR{1}(iit,:) = pr_mim;
        PR{2}(iit,:) = pr_mic;
        PR{3}(iit,:) = pr_aCoh;
        PR{4}(iit,:) = pr_iCoh;
        PR{5}(iit,:) = pr_absgc;
        PR{6}(iit,:) = pr_posgc;
        PR{7}(iit,:) = pr_abstrgc;
        PR{8}(iit,:) = pr_postrgc;
        
    catch
        a = [a iit];
    end
end

for ii = 1:length(PR)
    PR{ii}(a,:) = [];
end

%%

figure
figone(8,24)
ipip = 3;
oo = 1; 
for icon = [3 4 1 2 5 6 7 8 ] %Coherence, iCOH, MIM, MIC, absgc,posgc, abstrgc, postrgc
    
    data1 = PR{icon}(:,ipip);
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
    
 
%     grid on
    oo=oo+1;
    
    
end


outname = [DIRFIG 'figure2A.eps'];
print(outname,'-depsc');
% %%
% outname = [DIRFIG 'figure2A'];
% export_fig(outname, ['-r' '150'], '-a2', '-transparent'); 
close all











