DIRDATA = './';
DIRFIG = './figures/mimsim_ana/mim_sim2/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

% name = {...
%     'ip1_default';...
%     'ip2_int2';...
%     'ip2_int3';...
%     'ip2_int4';...
%     'ip2_int5';...
%     'ip3_iReg2';...
%     'ip4_snr005';...
%     'ip5_iss0';...
%     'ip5_iss025';...
%     'ip5_iss075';...
%     'ip5_iss1';...
%     'ip6_lag1';...
%     'ip7_dics';...
%     'ip7_eloreta';...
%     'ip8_hemisym1'};
% 
% 
% for inames = 1

%     clearvars -except mdefault DIRDATA DIRFIG name inames 
    
%     load([DIRDATA 'mim_sim1_auc_' name{inames} '.mat']) %auc = 100 iterations x 18 pipelines x 3 measures 

    load('mrr_mim2_default.mat')
    figure
    figone(30,50)
    
    if inames==1
        labs = {'mic','mim','mean icoh'};
    else
        labs = {'mic','mim','mean icoh','mean abscoh'};
    end

    for imim = 1

        subplot(2,1,imim)

        data = mrr(:,:,imim);
        data(data<0) = 0;
%         data(:,6,:)=[];

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
            xticklabels = {'1PC','2PCs', '3PCs','4PCs','5PCs',...
            '99%', '90%','sumVox','baseline'};
    
        else 
             xticklabels = {'1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%','baseline'};
        end
            
        
        set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
        ylabel(['MRR' labs(imim)])
        ylim([0 1])
       
        if imim==1 
            title(['SNR = 0.5, 1 interaction'])
        end

    end
    
%     if inames == 1 
%         %rearrange pipelines that they match reduced version 
%         mdefault(:,[1:8 10:11])=[];
%         mdefault1 = mdefault(:,2:end); 
%         mdefault1 = cat(2,mdefault1,mdefault(:,1)); 
%         mdefault=mdefault1; 
%         clear mdefault1
%     end
%     
%     outname = [DIRFIG name{inames} 'mrr'];
%     saveas(gcf,outname, 'png')
%     close all
% end 
