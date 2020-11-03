function fp_lfpmegmeg_ana 

DIRIN = '~/Dropbox/Franziska/Data_MEG_Project/lfpmegmeg/';
DIRFIG = '~/Dropbox/Franziska/Data_MEG_Project/figures/lfpmegmeg/';
if ~exist(DIRFIG);mkdir(DIRFIG); end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

for id = 1: numel(patientID)
    
    %% imCoh 
    clearvars -except DIRIN patientID id nit DIRFIG mim mic mic_shuf mim_shuf
    load(['MIM_GC_sub' patientID{id} '.mat'])
    nit = size(MIC_SHUF,1);
        
    
    for ii = 1:size(MIC_TRUE,1)
        for jj = 1: size(MIC_TRUE,2)
            for ifq = 1: size(MIC_TRUE,3)
                
                p_mic_permutation(ii,jj,ifq) = sum(MIC_TRUE(ii,jj,ifq)<MIC_SHUF(:,ii,jj,ifq))/nit;
                p_mim_permutation(ii,jj,ifq) = sum(MIM_TRUE(ii,jj,ifq)<MIM_SHUF(:,ii,jj,ifq))/nit;
                
            end
        end
    end
    
    mic(id,:,:,:) = MIC_TRUE;
    mim(id,:,:,:) = MIM_TRUE; 
    mim_shuf(id,:,:,:,:) = MIM_SHUF(1:100,:,:,:); 
    mic_shuf(id,:,:,:,:) = MIC_SHUF(1:100,:,:,:); 
    
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
    
    %%
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
%     caxis([4 15])
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
    
    %% GC
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
%     
    
end 
%%
[nsubs,nroi,~,nfreq] = size(granger);

%% true
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
%%
fbands = [1 4; 4 7; 7 15; 15 30; 30 45]; %theta, alpha, beta, low gamma, high gamma

for ibands = 1: size(fbands,1)
    figure
    subplot(1,2,1)
    imagesc(squeeze(sum(mask1(:,:,fbands(ibands,1):fbands(ibands,2)),3)))
    title('MIM')
    
    subplot(1,2,2)
    imagesc(squeeze(sum(mask2(:,:,fbands(ibands,1):fbands(ibands,2)),3)))
    title('MIC')
end

%%

