DIRDATA = './mim_sim4/';
DIRFIG = './figures/mimsim_ana/mim_sim4/';
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
[cb] = cbrewer2('spectral', 11);
cb1 = cbrewer2('Set1',9);

np = 6;

%%
for iname = 1%:numel(name)
    
    clearvars -except mdefault DIRDATA DIRFIG name iname labs cb cb1 np
    load([DIRDATA 'mrr_mim4_' name{iname} '.mat']);
    
    for imim = 1:4
        
        for imsr = 1:2
            if imsr==1
                if iname == 1
                    data = cat(2,mrr(:,np+7:end,imim),mrr(:,np+3:np+6,imim)); %zs0 + sumVox + baseline + 2 corrected pips
                else
                    data = mrr(:,:,imim);
                end
            elseif imsr == 2
                if iname==1
                    data = cat(2,pr(:,np+7:end,imim),pr(:,np+3:np+6,imim));
                else
                    data = pr(:,:,imim);
                end
            end
            
            if (imim == 1 || imim ==2)
                npips = size(data,2);
            elseif iname ==1 && imim>2
                data(:,np+3:np+4,:)=[];%no baseline, no sumVox, but corrected % pipelines
                npips = size(data,2);
            else
                npips = size(data,2)-1; %no baseline
            end
            %%
            figure
            figone(15,50)
            if iname == 1 && imim<3
                pips = [1:np+2 np+5 np+6 np+3 np+4];
            else
                pips = 1:npips;
            end
            count = 1;
            
            for ipip = pips
                
                data1 = squeeze(data(:,ipip));
                subplot(1,npips,count)
   
                if ipip<=np
                    cl = cb(np+1-ipip,:);
                elseif ipip <np+3
                    cl = cb(ipip+2,:);
                else
                    if iname ==1
                        if imim <3
                            if ipip == np+3
                                cl = cb1(8,:);
                            elseif ipip==np+4
                                cl = cb1(end,:);
                            elseif ipip > np+4
                                cl = cb(ipip-2,:);
                            end
                        else 
                            if ipip >np+2
                               cl = cb(ipip,:);
                            end
                        end 
                    else
                        cl = cb1(end,:);
                    end
                end
                [h, u] = fp_raincloud_plot(data1, cl, 1,0.2, 'ks');
                view([-90 -90]);
                set(gca, 'Xdir', 'reverse');
                set(gca, 'XLim', [0 1]);

                
                %                 grid on
                
                if iname == 1
                    if imim <3
                        xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','99%',...
                            '90%','sumVox','baseline','99% corrected','90% corrected'};
                    else
                        xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','99%',...
                            '90%','99% corrected','90% corrected'};
                    end
                else
                    xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','99%', '90%','baseline'};
                end
                
                title(xtitles(ipip))
                if count ~= 1
                    set(gca,'xtick',[])
                end
                
                if count == 1
                    if imsr==1
                        xlabel([labs{imim} ' MRR'])
                    else
                        xlabel([labs{imim} ' PR'])
                    end
                end
                
                %                 xlim([0 103])
                ylim([-0.75 2])
                count = count+1;
            end
            %%
            if imsr ==1
                outname = [DIRFIG name{iname} '_' labs{imim} '_mrr'];
            else
                outname = [DIRFIG name{iname} '_' labs{imim} '_pr'];
            end
            %saveas(gcf,outname, 'png')
            %close all
            
            
            %% ZS pipeline
            if iname ==1
                if imsr==1
                    data = mrr(:,1:np+2,imim); %zs pipelines
                elseif imsr == 2
                    data = pr(:,1:np+2,imim);
                end
                npips = size(data,2);
                
                figure
                figone(15,50)
                
                for ipip = 1:npips
                    
                    data1 = squeeze(data(:,ipip));
                    
                    subplot(1,npips,ipip)
                    if ipip<=np
                        cl = cb(np+1-ipip,:);
                    elseif ipip <np+3
                        cl = cb(ipip+2,:);
                    elseif ipip > np+4
                        cl = cb(ipip-2,:);
                    elseif ipip==np+4
                        cl = cb1(end,:);
                    else
                        cl = cb1(2,:);
                    end
                    [h, u] = fp_raincloud_plot(data1, cl, 1,0.2, 'ks');
                    view([-90 -90]);
                    set(gca, 'Xdir', 'reverse');
                    set(gca, 'XLim', [0 1]);
                    
                    xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','99%', '90%'};
                    title(xtitles(ipip))
                    
                    if ipip ~= 1
                        set(gca,'xtick',[])
                    end
                    
                    if ipip ==1 
                        if imsr==1
                            xlabel([labs{imim} ' MRR'])
                        else
                            xlabel([labs{imim} ' PR'])
                        end
                    end
                    
                    ylim([-0.75 2])
                    
                end
                
                if imsr ==1
                    outname = [DIRFIG 'ip9_ZS_' labs{imim} '_mrr'];
                else
                    outname = [DIRFIG 'ip9_ZS_' labs{imim} '_pr'];
                end
                
                %saveas(gcf,outname, 'png')
                %close all
            end
            %%
            
        end
        
    end
    
end


