function fp_hist_mrr_ICBEM

DIRDATA = './mim_sim4/';
DIRFIG = './figures/mimsim_ana/mim_sim4/ICBEM/';
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
    'ip7_champ_reg'};

labs = {'MIM','MIC','Mean abscoh','mean icoh','absGC','posGC','posGCw'};
%%
iname = 1;


clearvars -except iname name DIRDATA DIRFIG labs

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
    end
end


np = 6;

its = [1:nit];
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
        MRR{7}(iit,:) = mrr_posgc_w;
        
        PR{1}(iit,:) = pr_mim;
        PR{2}(iit,:) = pr_mic;
        PR{3}(iit,:) = pr_aCoh;
        PR{4}(iit,:) = pr_iCoh;
        PR{5}(iit,:) = pr_absgc;
        PR{6}(iit,:) = pr_posgc;
        PR{7}(iit,:) = pr_posgc_w;
        
        EM3{1}(iit,:) = em3_mim;
        EM3{2}(iit,:) = em3_mic;
        EM3{3}(iit,:) = em3_aCoh;
        EM3{4}(iit,:) = em3_iCoh;
        EM3{5}(iit,:) = em3_absgc;
        EM3{6}(iit,:) = em3_posgc;
        EM3{7}(iit,:) = em3_posgc_w;
        
        
        
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

for im = 3 %measures: MRR, PR, EM3
    
    for icon = 1:length(MRR) %MIM, MIC, aCoh, iCoh, absgc,posgc,posgc_w
        
        
        %%
        figure
        figone(15,50)
        
        if icon < 3 %MIC, MIM
            pips = [1:8 10 21 22 9]; %fixed, sumvox, summim, central, baseline
            npips = 12;
        else %GC
            pips = [1:8 21 22 9];
            npips = 11;
        end

        ipip1=1;
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
            
            subplot(1,npips,ipip1)
            
            
            [h, u] = fp_raincloud_plot(data1, cl, 1,0.2, 'ks');
            view([-90 -90]);
            set(gca, 'Xdir', 'reverse');
            set(gca, 'XLim', [0 1]);
            
            if icon < 3    
                xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                    '99%','MIM+sum','Sum+MIM','Central','TrueVox'};
            else
                xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                    '99%','Sum+MIM','Central','TrueVox'};
            end
            
            title(xtitles(ipip1))
            if ipip~=1
                set(gca,'xtick',[])
            end
            
            xlabel([labs{icon} ' ' imlab1])
            ipip1=ipip1+1;
            
        end
        
        ylim([-0.75 2])
        
        outname = [DIRFIG name{iname} '_' imlab '_' labs{icon}];
        saveas(gcf,outname, 'png')
        close all
        
        
    end
 
    
end












