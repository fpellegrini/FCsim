DIRDATA = './mimsim_ana/';
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

        subplot(numel(labs),1,imim)

        data = auc(:,:,imim);

        if inames == 1 || imim == 1 || imim ==2 
            npips = size(data,2); 
        else 
            npips = size(data,2)-1;
        end
        
        for ipip = 1:npips

            data1 = squeeze(data(:,ipip));
            md = mean(data1);
            if inames == 1 && (ipip == 9 || ipip >=12)
                mdefault(imim,ipip) = md;
            end
            a = ipip*3;

            ra = randi([-10 10],100,1)./20;
            scatter(ra+a,data1,'b')

            xlim([0 size(data,2)*10/3])
            if inames ~=1 && imim < 4 && ~(ipip==8 && imim==3)
                hold on
                plot([a-0.6 a+0.6],repmat(mdefault(imim,ipip),2,1),'g','linewidth',3)
            end
            hold on 
            if inames ==1 
                plot([a-0.6 a+0.6],repmat(md,2,1),'g','linewidth',3)
            else 
                plot([a-0.6 a+0.6],repmat(md,2,1),'r','linewidth',3)
            end
            
            grid on
            hold on
        end


        xTicks = (1:npips).*3;
        if inames == 1
            xticklabels = {'1zs','2zs', '3zs','4zs','5zs',...
            '99% zs','90% zs','sumVox','baseline','99%corr','90%corr', ...
            '1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%'};
    
        else 
             xticklabels = {'1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%','baseline'};
        end
            
        
        set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
        ylabel(['auc ' labs{imim}])
       
        if imim==1 
            title(name{inames}(5:end))
        end

    end
    
    if inames == 1 
        %rearrange pipelines that they match reduced version 
        mdefault(:,[1:8 10:11])=[];
        mdefault1 = mdefault(:,2:end); 
        mdefault1 = cat(2,mdefault1,mdefault(:,1)); 
        mdefault=mdefault1; 
        clear mdefault1
    end
    
    outname = [DIRFIG name{inames}];
    saveas(gcf,outname, 'png')
    close all
end 
