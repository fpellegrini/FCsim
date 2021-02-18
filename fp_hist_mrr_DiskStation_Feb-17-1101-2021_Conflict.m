DIRDATA = './mim_sim2/';
DIRFIG = './figures/mimsim_ana/mim_sim2/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

name = {...
    'ip1_default';...
    'ip2_int2';...
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
    'ip7_dics';...
    'ip7_eloreta_reg';...
    'ip8_hemisym1'};

labs = {'MIC','MIM','Mean icoh','mean abscoh'};
[cb] = cbrewer2('qual', 'Set1', 50, 'pchip');


%%
for iname = 1:5%numel(name)
    
    clearvars -except mdefault DIRDATA DIRFIG name iname labs
    load([DIRDATA 'mrr_mim2_' name{iname} '.mat']);
    
    for imim = 1:4
        
        for imsr = 1:2
            if imsr==1
                if iname == 1
                    data = cat(2,mrr(:,12:end,imim),mrr(:,8:11,imim)); %zs0 + sumVox + baseline + 2 corrected pips
                else
                    data = mrr(:,:,imim);
                end
            elseif imsr == 2
                if iname==1
                    data = cat(2,pr(:,12:end,imim),pr(:,8:11,imim));
                else
                    data = pr(:,:,imim);
                end
            end
            
            if (imim == 1 || imim ==2)
                npips = size(data,2);
            elseif iname ==1 && imim>2
                data(:,8:9,:)=[];%no baseline, no sumVox, but corrected % pipelines 
                npips = size(data,2);
            else
                npips = size(data,2)-1; %no baseline
            end
            
            figure
            figone(15,50)
            
            for ipip = 1:npips
                
                data1 = squeeze(data(:,ipip));
                subplot(1,npips,ipip)
                
%                 xbins =linspace(0.05,0.95,10);
%                 [counts,bins] = hist(data1,xbins);
                %  barh(bins,counts)
                
                [h, u] = raincloud_plot(data1, cb(3,:), 1,0.2, 'ks');
                view([-90 -90]);
                set(gca, 'Xdir', 'reverse');
                set(gca, 'XLim', [0 1]);

                grid on
                
                if iname == 1
                    if imim <3
                        xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','99%',...
                            '90%','sumVox','baseline','99% corrected','90% corrected'};
                    else
                        xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','99%',...
                            '90%','99% corrected','90% corrected'};
                    end
                else
                    xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','99%', '90%','baseline'};
                end
                
                title(xtitles(ipip))
                
                if imsr==1
                    ylabel([labs{imim} ' MRR'])
                else
                    ylabel([labs{imim} ' PR'])
                end
                
%                 xlim([0 103])
%                 ylim([0 1])
                
            end
            
            if imsr ==1
                outname = [DIRFIG name{iname} '_' labs{imim} '_mrr'];
            else
                outname = [DIRFIG name{iname} '_' labs{imim} '_pr'];
            end
            saveas(gcf,outname, 'png')
            close all
            
            
            %% ZS pipeline
            if iname ==1
                if imsr==1
                    data = mrr(:,1:7,imim); %zs pipelines 
                elseif imsr == 2
                    data = pr(:,1:7,imim);
                end
                npips = size(data,2);
                
                figure
                figone(15,50)
                
                for ipip = 1:npips
                    
                    data1 = squeeze(data(:,ipip));
                    xbins =linspace(0.05,0.95,10);
                    [counts,bins] = hist(data1,xbins);
                    subplot(1,npips,ipip)
                    barh(bins,counts)
                    grid on
                    
                    xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','99%', '90%'};
                    title(xtitles(ipip))
                    
                    if imsr==1
                        ylabel([labs{imim} ' MRR'])
                    else
                        ylabel([labs{imim} ' PR'])
                    end
                    
                    xlim([0 103])
                    ylim([0 1])
                    
                end
                
                if imsr ==1
                    outname = [DIRFIG 'ip9_ZS_' labs{imim} '_mrr'];
                else
                    outname = [DIRFIG 'ip9_ZS_' labs{imim} '_pr'];
                end
                
                saveas(gcf,outname, 'png')
                close all
            end
            %%
            
        end
        
    end
    
end


