function fp_hist_mrr_snr

DIRDATA = './mim_sim4/';
DIRFIG = './figures/mimsim_ana/mim_sim4/TRR/';
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

labs = {'MIM','MIC','Mean abscoh','mean icoh','absGC','posGC','posGCw'};

im = 2; %measures: MRR, PR, EM3    
icon = 1; %:length(MRR) %MIM, MIC, aCoh, iCoh, absgc,posgc,posgc_w
ipip = 3;
%%
o=1;
figure
figone(15,23)
for iname = [20 7 1 8]
    
    clearvars -except iname name DIRDATA DIRFIG labs o im icon ipip
    
    %default paramenters
    nit = 100;
    iInt = 2;
    iReg=1;
    isnr=0.7;
    iss = 0.5;
    ilag=2;
    ifilt='l';
    
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
                dimred = 'pr';
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
            case 20
                dimred = 'p';
                isnr = 0.3;
        end
    end
    
    
    np = 6;
    
    its = [1:nit];
    a=[];
    
    %%
    for iit= its
        
        try
            if iname == 18 || iname == 17 || iname == 20 || iname == 7
                inname = sprintf('mrr_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d_%s'...
                    ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit,dimred);
            else
                inname = sprintf('mrr_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
                    ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
            end
            
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
            
            EM3{1}(iit,:) = em3_mim;
            EM3{2}(iit,:) = em3_mic;
            EM3{3}(iit,:) = em3_aCoh;
            EM3{4}(iit,:) = em3_iCoh;
            EM3{5}(iit,:) = em3_absgc;
            EM3{6}(iit,:) = em3_posgc;
            
            
            
        catch
            a = [a iit];
        end
    end
    
    for ii = 1:length(MRR)
        MRR{ii}(a,:) = [];
        PR{ii}(a,:) = [];
        EM3{ii}(a,:) = [];
    end
        
%%
    
    switch im
        case 1
            data1 = MRR{icon}(:,ipip);
            imlab = 'MRR';
            imlab1 = 'MRR';
        case 2
            data1 = PR{icon}(:,ipip);
            imlab = 'PR';
            imlab1 = 'PR';
        case 3
            data1 = EM3{icon}(:,ipip);
            imlab = 'EM';
            imlab1 = '1-EM';
    end
    
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
    
    if im == 2
        [h, u] = fp_raincloud_plot_a(data1, cl, 1,0.2, 'ks');
    else
        [h, u] = fp_raincloud_plot(data1, cl, 1,0.2, 'ks');
    end
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    set(gca, 'XLim', [0 1]);
    
    titles = {'snr 0.3', 'snr 0.5','snr 0.7','snr 0.9'};
    title(titles{o})
    ylim([-0.75 2])
    if o==1
        xlabel([labs{icon} ' ' imlab1])
    end
    grid on
    set(gca,'ytick',[])
%     
%     if o~=1
%         set(gca,'xtick',[])
%     end
    
    o=o+1;
    
end








outname = [DIRFIG imlab '_' labs{icon} '_snr_3pcs'];
saveas(gcf,outname, 'png')
close all












