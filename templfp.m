
DIRIN = '~/Dropbox/Franziska/Data_MEG_Project/lfpmegmeg/5fixed_noZS/';
DIRFIG = '~/Dropbox/Franziska/Data_MEG_Project/figures/lfpmegmeg/5fixed_noZS/';
if ~exist(DIRFIG);mkdir(DIRFIG); end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
scores = db_twstrs();
t_score = nan(numel(patientID),1);

%%
for id = 1: numel(patientID)
    
    % imCoh (mim/mic) 
    clearvars -except DIRIN patientID id nit DIRFIG mim mic mic_shuf mim_shuf t_score scores
    
     if ~isnan(scores{1,id}{1,2})
         t_score(id) = scores{1,id}{1,2};
     end
    
    
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
    
    % GC
    
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
end
%%
fbands = [1 4; 4 7; 7 15; 15 30; 30 45]; %theta, alpha, beta, low gamma, high gamma
fnames = {'theta', 'alpha', 'beta', 'low gamma', 'high gamma'};

for ibands = 1: size(fbands,1)
   
    clear mimf_s micf_s 
    gf(:,:,:,ibands) = mean(granger(:,:,:,fbands(ibands,1):fbands(ibands,2)),4);
    mimf(:,:,:,ibands) = mean(mim(:,:,:,fbands(ibands,1):fbands(ibands,2)),4);
    micf(:,:,:,ibands) = mean(mic(:,:,:,fbands(ibands,1):fbands(ibands,2)),4);
    
    mimf_s = mean(mim_shuf(:,:,:,:,fbands(ibands,1):fbands(ibands,2)),5);
    micf_s = mean(mic_shuf(:,:,:,:,fbands(ibands,1):fbands(ibands,2)),5);
    
    [p_gc(:,:,ibands),~,~,~] = fp_get_signrank_results_gc(gf(:,:,:,ibands),0.001);
    [~, mask_gc(:,:,ibands)] = fdr( p_gc(:,:,ibands), 0.05);

    [p_mim(:,:,ibands),~,~] = fp_get_signrank_results_megmeg(mimf(:,:,:,ibands),mimf_s,0.001);
    [p1, mask_mim(:,:,ibands)] = fdr(p_mim(:,:,ibands), 0.05);

    [p_mic(:,:,ibands),~,~] = fp_get_signrank_results_megmeg(micf(:,:,:,ibands),micf_s,0.001);
    [p2, mask_mic(:,:,ibands)] = fdr(p_mic(:,:,ibands), 0.05);
    
end
% outname= [DIRFIG 'MIM_group_band'];
% saveas(gcf,outname, 'png')
% close all
%%
mask_mim=double(mask_mim);
mask_mim(mask_mim==0) = 0.3;
mask_mic=double(mask_mic);
mask_mic(mask_mic==0) = 0.3;
mask_gc=double(mask_gc);
mask_gc(mask_gc==0) = 0.3;

figone(30,60)

for ibands = 1:5
    
    subplot(3,5,ibands)
    imagesc(squeeze(-log10(p_mim(:,:,ibands))),'AlphaData',mask_mim(:,:,ibands))
    title(['-log10p MIM ', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')  
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45) 
    caxis([0 4])
    
    subplot(3,5,ibands+5)
    imagesc(squeeze(-log10(p_mic(:,:,ibands))),'AlphaData',mask_mic(:,:,ibands))
    title(['-log10p MIC ', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')  
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45) 
    caxis([0 4])
    
    subplot(3,5,ibands+10)
    imagesc(squeeze(-log10(p_gc(:,:,ibands))),'AlphaData',mask_gc(:,:,ibands))
    title(['-log10p GC ', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')  
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45) 
    caxis([0 4])
    if ibands == 5
        colorbar
    end
    
end

outname= [DIRFIG 'Group_analysis'];
saveas(gcf,outname, 'png')

%%

inds = find(isnan(t_score));
t_score(inds)=[];
mimf(inds,:,:,:)=[];
micf(inds,:,:,:)=[];
gf(inds,:,:,:)=[];
nroi=size(gf,2);


for ibands = 1:5
    for iroi = 1:nroi
        for jroi = 1:nroi
            [rho_mimf(iroi,jroi,ibands), p_tmim(iroi,jroi, ibands)] = corr(t_score, mimf(:,iroi,jroi,ibands), 'Type','Spearman','tail','right');
            [rho_micf(iroi,jroi,ibands), p_tmic(iroi,jroi, ibands)] = corr(t_score, micf(:,iroi,jroi,ibands), 'Type','Spearman','tail','right');
            [rho_gf(iroi,jroi,ibands), p_tg(iroi,jroi, ibands)] = corr(t_score, gf(:,iroi,jroi,ibands), 'Type','Spearman');
            
        end
    end
    [~, mask_tmim(:,:,ibands)] = fdr(p_tmim(:,:,ibands),0.05);
    [~, mask_tmic(:,:,ibands)] = fdr(p_tmic(:,:,ibands),0.05);
    [~, mask_tgc(:,:,ibands)] = fdr(p_tg(:,:,ibands),0.05);
end

%%
mask_tmim=double(mask_tmim);
mask_tmim(mask_tmim==0) = 0.3;
mask_tmic=double(mask_tmic);
mask_tmic(mask_tmic==0) = 0.3;
mask_tgc=double(mask_tgc);
mask_tgc(mask_tgc==0) = 0.3;

figone(30,60)

for ibands = 1:5
    
    subplot(3,5,ibands)
    imagesc(squeeze(-log10(p_tmim(:,:,ibands))),'AlphaData',mask_tmim(:,:,ibands))
    title(['-log10p MIM corr ', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')  
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45) 
    caxis([0 4])
    
    subplot(3,5,ibands+5)
    imagesc(squeeze(-log10(p_tmic(:,:,ibands))),'AlphaData',mask_tmic(:,:,ibands))
    title(['-log10p MIC corr', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')  
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45) 
    caxis([0 4])
    
    subplot(3,5,ibands+10)
    imagesc(squeeze(-log10(p_tg(:,:,ibands))),'AlphaData',mask_tgc(:,:,ibands))
    title(['-log10p GC corr', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')  
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45) 
    caxis([0 4])
    if ibands == 5
        colorbar
    end
    
end

outname= [DIRFIG 'Group_analysis_TWSTRS'];
saveas(gcf,outname, 'png')
    