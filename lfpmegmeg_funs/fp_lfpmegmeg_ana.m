function fp_lfpmegmeg_ana 

DIRIN = '~/Dropbox/Franziska/Data_MEG_Project/lfpmegmeg/5fixed_noZS/';
DIRFIG = '~/Dropbox/Franziska/Data_MEG_Project/figures/lfpmegmeg/5fixed_noZS/';
if ~exist(DIRFIG);mkdir(DIRFIG); end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

%%
for id = 1: numel(patientID)
    
    %% imCoh (mim/mic) 
    clearvars -except DIRIN patientID id nit DIRFIG mim mic mic_shuf mim_shuf
    load(['MIM_GC_sub' patientID{id} '.mat'])
    nit = size(MIC_SHUF,1);
    
    % p values of permutation testing per sub, roi-roi combi and freq 
    for ii = 1:size(MIC_TRUE,1)
        for jj = 1: size(MIC_TRUE,2)
            for ifq = 1: size(MIC_TRUE,3)
                
                p_mic_permutation(ii,jj,ifq) = sum(MIC_TRUE(ii,jj,ifq)<MIC_SHUF(:,ii,jj,ifq))/nit;
                p_mim_permutation(ii,jj,ifq) = sum(MIM_TRUE(ii,jj,ifq)<MIM_SHUF(:,ii,jj,ifq))/nit;
                
            end
        end
    end
    
    %collect mim and mic values of every sub 
    mic(id,:,:,:) = MIC_TRUE;
    mim(id,:,:,:) = MIM_TRUE; 
    mim_shuf(id,:,:,:,:) = MIM_SHUF(1:100,:,:,:); 
    mic_shuf(id,:,:,:,:) = MIC_SHUF(1:100,:,:,:); 
    
    %% GC
    
    % indices of DIFFGC_TRUE
    ninds=1;
    for ipcs = 1:10
        for jpcs = ipcs+1:10
            inds(ninds,:) = [ipcs,jpcs];
            ninds= ninds+1;
        end
    end
    
    granger1 = zeros(10,10,size(DIFFGC_TRUE,2)); 
    for iind = 1: size(inds,1)
        granger1(inds(iind,1),inds(iind,2),:) = DIFFGC_TRUE(iind,:);
    end
    
    granger(id,:,:,:) = granger1; 
    
    %% subject-wise plots of true mim/mic
%     figone(40,40)
%     subplot(2,2,1)
%     imagesc(squeeze(sum(MIC_TRUE,3)))
%     title('MIC, sum across freqs')
%     xlabel('Regions')
%     ylabel('Regions')
%         xticks=1:10;
%     xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
%         'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
%     set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',xticks,'YTickLabel',xticklabels)
%     xtickangle(45)
%     colorbar
%     caxis([4 22])
% 
%     subplot(2,2,2)
%     imagesc(squeeze(MIC_TRUE(:,7,:)))
%     title('MIC at cerebellum')
%     xlabel('Frequency in Hz')
%     ylabel('Regions')  
%     xticks = 1:5:45;
%     xticklabels=0:10:90;
%     yticks=1:10;
%     yticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
%         'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};    
%     set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels)
%     colorbar
%     caxis([0 0.4])
%     
%     
%     
%     subplot(2,2,3)
%     a = squeeze(-log10(p_mic_permutation(:,:,13)));
%     a(a<1.3)=0;
%     imagesc(a)
%     title('-log10(pval) MIC at 26 Hz')
%     xlabel('Regions')
%     ylabel('Regions')  
%     xticks=1:10;
%     xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
%         'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
%     set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',xticks,'YTickLabel',xticklabels)
%     xtickangle(45)
%     colorbar
%     caxis([0 2])
%     
%     
%     subplot(2,2,4)
%     b = squeeze(-log10(p_mic_permutation(:,7,:)));
%     b(b<1.3)=0;
%     imagesc(b)
%     title('-log10(pval) MIC at cerebellum')
%     xlabel('Frequency in Hz')
%     ylabel('Regions')  
%     xticks = 1:5:45;
%     xticklabels=0:10:90;
%     yticks=1:10;
%     yticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
%         'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};    
%     set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels)
%     colorbar    
%     caxis([0 2])
%     
%     outname= [DIRFIG 'imCoh_sub' patientID{id}];
%     saveas(gcf,outname, 'png')
%     close all
%     
%     % subject-wise plots of true GC
%     roi = 7;
%     
%     figone(15,40) 
%     subplot(1,2,1)
%     imagesc(sum(granger1,3))
%     title('DIFFGC, sum across freqs')
%     xlabel('Regions')
%     ylabel('Regions')
%     xticks=1:10;
%     xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
%         'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
%     set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',xticks,'YTickLabel',xticklabels)
%     xtickangle(45)
%     colorbar
%     caxis([-4 4])
%     
%     
%     subplot(1,2,2)
%     cer = squeeze(granger1(:,roi,:));
%     cer(roi:end,:) = squeeze(granger1(roi,roi:end,:));
%     imagesc(cer)
%     title('DIFFGC at cerebellum')
%     xlabel('Frequency in Hz')
%     ylabel('Regions')  
%     xticks = 1:5:45;
%     xticklabels=0:10:90;
%     yticks=1:10;
%     yticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
%         'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};    
%     set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels)
%     colorbar    
%     caxis([-0.25 0.25])
%     
%     outname= [DIRFIG 'GC_sub' patientID{id}];
%     saveas(gcf,outname, 'png')
%     close all
    
    
end 

[nsubs,nroi,~,nfreq] = size(granger);

%% between subjects for 5 freq bands 
fprintf('Testing...\n')

[p_gc,~,~,~] = fp_get_signrank_results_gc(granger,0.001);
[p, mask] = fdr( p_gc, 0.05);
sum(mask(:))

[p_mim,~,~] = fp_get_signrank_results_megmeg(mim,mim_shuf,0.001);
[p1, mask1] = fdr(p_mim, 0.05);
sum(mask1(:))

[p_mic,~,~] = fp_get_signrank_results_megmeg(mic,mic_shuf,0.001);
[p2, mask2] = fdr(p_mic, 0.05);
sum(mask2(:))
%
fbands = [1 4; 4 7; 7 15; 15 30; 30 45]; %theta, alpha, beta, low gamma, high gamma
fnames = {'theta', 'alpha', 'beta', 'low gamma', 'high gamma'};

figone(30,60)
for ibands = 1: size(fbands,1)
    
    subplot(2,5,ibands)
    imagesc(squeeze(sum(mask1(:,:,fbands(ibands,1):fbands(ibands,2)),3))./length(fbands(ibands,1):fbands(ibands,2)))
    title(['MIM ', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')  
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45) 
    caxis([0 1])
    
    subplot(2,5,ibands+5)
    imagesc(squeeze(sum(mask2(:,:,fbands(ibands,1):fbands(ibands,2)),3))./length(fbands(ibands,1):fbands(ibands,2)))
    title(['MIC ', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')  
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45)
    caxis([0 1])
    
end
outname= [DIRFIG 'MIM_group_band'];
saveas(gcf,outname, 'png')
close all

%% within subjects for 5 freq bands 

for iband = 1: size(fbands,1)   
    mim_bands(:,:,:,iband) = squeeze(sum(mim(:,:,:,fbands(ibands,1):fbands(ibands,2)),4));
    mic_bands(:,:,:,iband) = squeeze(sum(mic(:,:,:,fbands(ibands,1):fbands(ibands,2)),4));
    
    mim_shuf_bands(:,:,:,:,ibands) = squeeze(sum(mim_shuf(:,:,:,:,fbands(ibands,1):fbands(ibands,2)),5));
    mic_shuf_bands(:,:,:,:,ibands) = squeeze(sum(mic_shuf(:,:,:,:,fbands(ibands,1):fbands(ibands,2)),5));
end

for isub = 1: size(mim,1)
    for iroi = 1:size(mim,2)
        for jroi = 1:size(mim,3)
            for iband = 1: size(fbands,1)
                
                p_mim_bands(isub,iroi,jroi,iband) = sum(mim_shuf_bands(isub, :,iroi,jroi,iband)>mim_bands(isub,iroi,jroi,iband))/nit;
                p_mic_bands(isub,iroi,jroi,iband) = sum(mic_shuf_bands(isub, :,iroi,jroi,iband)>mic_bands(isub,iroi,jroi,iband))/nit;
            end
        end
    end
end

logp_mim = -log10(p_mim_bands); 
logp_mim(logp_mim<=1.3)=0;
logp_mim(isinf(logp_mim))=2;
logp_mic = -log10(p_mic_bands); 
logp_mic(logp_mic<=1.3)=0;
logp_mic(isinf(logp_mic))=2;


for iband = 1:size(fbands,1)
    figure 
    figone(30,60)
    for isub = 1:size(mim,1) 

        subplot(3,4,isub) 
        imagesc(squeeze(logp_mim(isub,:,:,iband)))
        title(['MIM ', fnames{iband}])
        xlabel('Regions')
        ylabel('Regions')  
        xticks=1:10;
        xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
            'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
        set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',xticks,'YTickLabel',xticklabels)
        xtickangle(45)
        colorbar
        caxis([0 2])
        
        
    end 
    
    outname= [DIRFIG 'MIM_withinsubs_band' fnames{iband}];
    saveas(gcf,outname, 'png')
    close all
end 

for iband = 1:size(fbands,1)
    figure 
    figone(30,60)
    for isub = 1:size(mim,1) 
        
        subplot(3,4,isub)
        imagesc(squeeze(logp_mic(isub,:,:,iband)))
        title(['MIC ', fnames{iband}])
        xlabel('Regions')
        ylabel('Regions')  
        xticks=1:10;
        xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
            'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
        set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',xticks,'YTickLabel',xticklabels)
        xtickangle(45)
        colorbar
        caxis([0 2])
    end 
    
    outname= [DIRFIG 'MIC_withinsubs_band' fnames{iband}];
    saveas(gcf,outname, 'png')
    close all
end 
