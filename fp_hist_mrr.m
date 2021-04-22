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

labs = {'MIC','MIM','Mean abscoh','mean icoh','absGC','posGC','posGCw'};
[cb] = cbrewer2('spectral', 11);
cb1 = cbrewer2('Set1',9);

%%

for iname = 1
    
    
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
            MRR{7}(iit,:) = mrr_posgc_w;
            
            PR{1}(iit,:) = pr_mim;
            PR{2}(iit,:) = pr_mic;
            PR{3}(iit,:) = pr_aCoh;
            PR{4}(iit,:) = pr_iCoh;
            PR{5}(iit,:) = pr_absgc;
            PR{6}(iit,:) = pr_posgc;
            PR{7}(iit,:) = pr_posgc_w;
            
            EM1{1}(iit,:) = em1_mim;
            EM1{2}(iit,:) = em1_mic;
            EM1{3}(iit,:) = em1_aCoh;
            EM1{4}(iit,:) = em1_iCoh;
            EM1{5}(iit,:) = em1_absgc;
            EM1{6}(iit,:) = em1_posgc;
            EM1{7}(iit,:) = em1_posgc_w;
             
            EM2{1}(iit,:) = em2_mim;
            EM2{2}(iit,:) = em2_mic;
            EM2{3}(iit,:) = em2_aCoh;
            EM2{4}(iit,:) = em2_iCoh;
            EM2{5}(iit,:) = em2_absgc;
            EM2{6}(iit,:) = em2_posgc;
            EM2{7}(iit,:) = em2_posgc_w;
            
                        
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
    
    for ii = 1:6
        MRR{ii}(a,:) = [];
        PR{ii}(a,:) = [];
    end
    
    %%
    
    for im = 1:5 %measures: MRR, PR, 3 EM measures 
        
        for icon = 1:length(MRR) %MIC, MIM, aCoh, iCoh, absgc,posgc,posgc_w 


            %%
            figure
            if iname == 15
                figone(15,25)
            else
                figone(15,50) 
            end

            if icon < 3 %MIC, MIM
                if iname == 1
                    pips = 1:12;
                    npips = 12;
                else 
                    pips = 1:9;
                    npips = 9;
                end
            elseif icon > 2 && icon < 5
                 if iname == 1
                    pips = [1:9 11:12];%exclude sumVox pip
                    npips = 11;
                else 
                    pips = 1:9;
                    npips = 9;
                 end
            else
                pips = 1:9;
                npips = 9;
            end

            if iname == 15 %champ
                pips = [1:3 8];
                npips = 4;
            end

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
                        data1 = EM1{icon}(:,ipip);
                        imlab = 'EM1';
                        imlab1 = '1-EM1';
                    case 4
                        data1 = EM2{icon}(:,ipip);
                        imlab = 'EM2';
                        imlab1 = '1-EM2';
                    case 5
                        data1 = EM3{icon}(:,ipip);
                        imlab = 'EM3';
                        imlab1 = '1-EM3';
                end 

                if ipip<=np
                    cl = [0.5 0.5 1];%cb(np+1-ipip,:);
                elseif ipip <np+3 || ipip == 11 || ipip == 12
                    cl = [1 0.5 0.5];%cb(ipip+2,:);
                elseif ipip == 10
                        cl = [0.5 1 0.5];
                elseif ipip == 9
                    cl = [0.2 0.2 0.2];
                end

                if iname == 15 && ipip == 8 
                    ipip1 = 4;
                elseif iname == 1 && icon > 2 && icon < 5 && ipip > 9
                    ipip1 = ipip-1;
                else 
                    ipip1 = ipip;
                end
                subplot(1,npips,ipip1)

                [h, u] = fp_raincloud_plot(data1, cl, 1,0.2, 'ks');
                view([-90 -90]);
                set(gca, 'Xdir', 'reverse');
                set(gca, 'XLim', [0 1]);
                
                if iname == 1
                    
                    xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                        '99%','baseline','sumVox','90% corrected','99% corrected'};


%                 elseif iname == 15 
%                     xtitles = {'1PC','2PCs', '3PCs','99%'};
                else 
                    xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%', '99%','baseline'};
                end

                title(xtitles(ipip))
                if ipip~=1
                    set(gca,'xtick',[])
                end

                xlabel([labs{icon} ' ' imlab1])

                if iname == 15 && ipip == 4
                    break
                end

            end

            %                 xlim([0 103])
            ylim([-0.75 2])

            outname = [DIRFIG name{iname} '_' imlab '_' labs{icon}];
            saveas(gcf,outname, 'png')
            close all


        end


        if iname == 1

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
                    
                    switch im
                        case 1
                            data1 = MRR{icon}(:,ipip+12);
                            imlab = 'MRR';
                            imlab1 = 'MRR';
                        case 2
                            data1 = PR{icon}(:,ipip+12);
                            imlab = 'PR';
                            imlab1 = 'PR';
                        case 3
                            data1 = EM1{icon}(:,ipip+12);
                            imlab = 'EM1';
                            imlab1 = '1-EM1';
                        case 4
                            data1 = EM2{icon}(:,ipip+12);
                            imlab = 'EM2';
                            imlab1 = '1-EM2';
                        case 5
                            data1 = EM3{icon}(:,ipip+12);
                            imlab = 'EM3';
                            imlab1 = '1-EM3';
                    end
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


                    xlabel([labs{icon} ' ' imlab1])

                end

                %                 xlim([0 103])
                ylim([-0.75 2])

                outname = [DIRFIG name{iname} '_' imlab '2_' labs{icon}];
                saveas(gcf,outname, 'png')
                close all


            end

        end
    end
    
   
    
    
end




    
    

