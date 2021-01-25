DIRDATA = './';
DIRFIG = './figures/mimsim_ana/';

name = {...
    'ip1_default';...
    'ip2_int2';...
    'ip2_int3';...
    'ip2_int4';...
    'ip2_int5';...
    'ip3_iReg2';...
    'ip4_snr005';...
    'ip5_iss0';...
    'ip5_iss025';...
    'ip5_iss075';...
    'ip5_iss1';...
    'ip6_lag1';...
    'ip7_dics';...
    'ip7_eloreta';...
    'ip8_hemisym1'};


for inames = 1:numel(name)
    
    clearvars -except mdefault DIRDATA DIRFIG name inames
    
    load([DIRDATA 'mim_sim1_auc_' name{inames} '.mat']) %auc = 100 iterations x 18 pipelines x 3 measures
    
    figure
    figone(30,50)
    
    if inames==1
        labs = {'mic','mim','mean icoh'};
    else
        labs = {'mic','mim','mean icoh','mean abscoh'};
    end
    
    for imim = 1:numel(labs)
        
        %         subplot(numel(labs),1,imim)
        
        data = auc(:,:,imim);
        
        npips = size(data,2);
        %         if inames == 1 || imim == 1 || imim ==2
        %             npips = size(data,2);
        %         else
        %             npips = size(data,2)-1;
        %         end
        
        for ipip = 1:npips
            
            data1 = squeeze(data(:,ipip));
            %             md = mean(data1);
            %             if inames == 1 && (ipip == 9 || ipip >=12)
            %                 mdefault(imim,ipip) = md;
            %             end
            
            if ~(imim >2 && ipip == 8)
                [counts,bins] = hist(data1,10);
                subplot(numel(labs),npips,ipip+(imim*npips)-npips)
                barh(bins,counts)
                
                xTicks = (1:npips).*3;
                if inames == 1
                    xticklabels = {'1zs','2zs', '3zs','4zs','5zs',...
                        '99% zs','90% zs','sumVox','baseline','99%corr','90%corr', ...
                        '1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%'};
                    
                else
                    xticklabels = {'1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%','baseline'};
                end
                
                if imim ==1
                    title(xticklabels(ipip))
                end
                xlabel(['auc ' labs{imim}])
                xlim([0 30])
            end
            
        end
        
    end
    
    
    outname = [DIRFIG name{inames} '_hist'];
    saveas(gcf,outname, 'png')
    close all
end
