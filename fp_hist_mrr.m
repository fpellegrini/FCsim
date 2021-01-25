DIRDATA = './';
DIRFIG = './figures/mimsim_ana/mim_sim2/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

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
    'ip8_hemisym1';...
    'ip_4_snr01'};
%
%
% for inames = 1
inames=13;
%     clearvars -except mdefault DIRDATA DIRFIG name inames

%     load([DIRDATA 'mim_sim1_auc_' name{inames} '.mat']) %auc = 100 iterations x 18 pipelines x 3 measures

load('mrr_mim2_dics.mat')
figure
figone(30,50)

if inames==1
    labs = {'MIC','MIM','Mean icoh'};
else
    labs = {'MIC','MIM','Mean icoh','mean abscoh'};
end
count = 1;
for imim = [1 3]
    
%     subplot(2,1,count)
    
    data = mrr(:,:,imim);
    data(data<0) = 0;
    %         data(:,6,:)=[];
    
    if  (imim == 1 || imim ==2)
        npips = size(data,2);
        npips1 = size(data,2);
    elseif inames==1 && imim ==3
        npips = size(data,2)-2;
        npips1 = size(data,2);
    elseif inames>1 && imim ==3 
        npips = size(data,2)-1;
        npips1 = size(data,2);
    else
        npips = size(data,2);
        npips1 = size(data,2);
    end
    
    for ipip = 1:npips
        
        data1 = squeeze(data(:,ipip));
        [counts,bins] = hist(data1,10);
        subplot(2,npips1,ipip+(count*npips1)-npips1)
        barh(bins,counts)
        grid on
        
        xTicks = (1:npips).*3;
        if inames == 1
            xticklabels = {'1PC','2PCs', '3PCs','4PCs','5PCs',...
                '99%', '90%','sumVox','baseline'};
            
        else
            xticklabels = {'1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%','baseline'};
        end

        title(xticklabels(ipip))

        ylabel([labs{imim} ' MRR'])
        xlim([0 100])
        ylim([0 1])
        
    end
    count = count +1; 
end

outname = [DIRFIG name{inames} 'mrr'];
saveas(gcf,outname, 'png')


