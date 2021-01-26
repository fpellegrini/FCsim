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
    'ip7_eloreta';...
    'ip8_hemisym1'};

labs = {'MIC','MIM','Mean icoh','mean abscoh'};



for iname = 1:numel(name)
    
    clearvars -except mdefault DIRDATA DIRFIG name iname labs
    load([DIRDATA 'mrr_mim2_' name{iname} '.mat']);
      
    %     count = 1;
    for imim = 1:3
        
        for imsr = 1:3
            if imsr==1
                data = mrr(:,:,imim);
            elseif imsr == 2
                data = pr(:,:,imim);
            end
            
            if (imim == 1 || imim ==2)
                npips = size(data,2);
                %                 npips1 = size(data,2);
            else
                npips = size(data,2)-1;
                %             elseif inames==1 && imim ==3
                %                 npips = size(data,2)-2;
                %                 npips1 = size(data,2);
                %             elseif inames>1 && imim ==3
                %                 npips = size(data,2)-1;
                %                 npips1 = size(data,2);
                %             else
                %                 npips = size(data,2);
                %                 npips1 = size(data,2);
            end
            
            figure
            figone(15,50)
            
            for ipip = 1:npips
                
                data1 = squeeze(data(:,ipip));
                xbins =linspace(0.05,0.95,10); %linspace(0,1,10);
                [counts,bins] = hist(data1,xbins);
                %                 subplot(2,npips1,ipip+(count*npips1)-npips1)
                subplot(1,npips,ipip)
                barh(bins,counts)
                grid on
                
                xTicks = (1:npips).*3;
                if iname == 1
                    xticklabels = {'1PC','2PCs', '3PCs','4PCs','5PCs','99%', '90%','sumVox','baseline'};
%                     xticklabels = {'1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%','baseline'};
                    
                else
                    xticklabels = {'1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%','baseline'};
                end
                
                title(xticklabels(ipip))
                
                if imsr==1
                    ylabel([labs{imim} ' MRR'])
                else
                    ylabel([labs{imim} ' PR'])
                end
                
                xlim([0 103])
                ylim([0 1])
                
            end
            %             count = count +1;
            
            if imsr ==1
                outname = [DIRFIG name{iname} '_' labs{imim} '_mrr'];
            else
                outname = [DIRFIG name{iname} '_' labs{imim} '_pr'];
            end
            saveas(gcf,outname, 'png')
            close all
            
        end
              
    end
    
end


