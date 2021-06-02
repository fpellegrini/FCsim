function fp_corrs_ICBEM

DIRDATA = './mim_sim4/';
DIRFIG = './figures/mimsim_ana/mim_sim4/ICBEM/';
if ~exist(DIRFIG); mkdir(DIRFIG); end

labs = {'MIC','MIM','absCOH','iCOH','npcs'};

%% only for default setting 

clearvars -except iname name DIRDATA DIRFIG labs

%default paramenters
nit = 100;
iInt = 2;
iReg=1;
isnr=0.7;
iss = 0.5;
ilag=2;
ifilt='l';

its = [1:nit];
a=[];

%%
for iit= its
    
    try
        inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
            ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
        
        load([DIRDATA inname '.mat'])
        
        CORR(iit,:,1) = corr_voxmim;
        CORR(iit,:,2) = corr_voxmic;
        CORR(iit,:,3) = corr_voxacoh;
        CORR(iit,:,4) = corr_voxicoh;
        CORR(iit,:,5) = corr_voxnpcs;
        
        for ipip = 1:8
            VAREX(iit,ipip,:) = to_save{ipip}.varex;
        end
        
    catch
        a = [a iit];
    end
end

CORR1 = CORR(:,1:8,:);
CORR1(a,:,:)=[];
VAREX(a,:,:)=[];

%% corrs 

figone(30,30)

for icon = 1:size(CORR1,3)-1 %MIC, MIM, aCoh, iCoh, npcs
    
    subplot(2,2,icon)
        
    p1 = bar(1:6,squeeze(mean(CORR1(:,1:6,icon),1)));    
    hold on
    p2 = bar(7:8,squeeze(mean(CORR1(:,7:8,icon),1)));  
    
    set(p1,'FaceColor',[0.8 0.7 0.6])
    set(p2,'FaceColor',[0.8 0.4 0.5])
    
    xticklabels = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
        '99%'};
    
    set(gca,'xtick',1:8,'xticklabels',xticklabels)
    xtickangle(45)
    title(['nvoxels per region x ' labs{icon}])
    ylabel('Correlation coefficient')
    ylim([-0.8 0.8])
    
    
end

outname = [DIRFIG 'corrs_default'];
saveas(gcf,outname, 'png')
close all

%% varex 

figone(15,30)
var1 = mean(VAREX,3);
for ii = 1:size(VAREX,2)
    subplot(1,8,ii)
    
    if ii <= 6 
        cl = [0.8 0.7 0.6];
    else 
        cl = [0.8 0.4 0.5];
    end
    
    [h, u] = fp_raincloud_plot(var1(:,ii), cl, 1,0.2, 'ks');
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    set(gca, 'XLim', [0 1]);
    title(xticklabels{ii})
    if ii ==1
        xlabel('Variance explained')
    else
        set(gca,'xtick',[])
    end
end

outname = [DIRFIG 'varex_default'];
saveas(gcf,outname, 'png')
close all










