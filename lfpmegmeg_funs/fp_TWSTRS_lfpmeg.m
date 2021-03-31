function fp_TWSTRS_lfpmeg

DIRDATA = './lfpmegmeg/';

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

scores = db_twstrs();
alpha = 0.05;

nsubs = numel(patientID);
%% collect twstrs scores, MIC and GC of all subjects 
o = 1;
for id = 1:nsubs
    clearvars -except DIRDATA patientID scores alpha nsubs o id t_score tCoh granger
    
    if ~isnan(scores{1,id}{1,2})
        fprintf('Working on subject %s \n',patientID{id})
        
        t_score(o) = scores{1,id}{1,2};
        load([DIRDATA 'MIM_GC_sub' patientID{id} '.mat'])
        
        tCoh(o,:,:,:) = MIC_TRUE; %nsubs x nrois x nrois x nfreqs
        %sCoh(o,:,:,:,:) = MIC_SHUF; %nsubs x nit x nrois x nrois x nfreqs
        
        % granger
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
        
        granger(o,:,:,:) = granger1;
    end
    o = o+1;
end

[nsubs, nroi, ~, nfreq] = size(tCoh);
granger(:,:,:,1) = []; %remove offset frequency

%% correlation with TWSTRS 

for ifreq = 1:nfreq
    for iroi = 1:nroi
        for jroi = 1:nroi
            [rho_c(iroi,jroi,ifreq), pval_c(iroi,jroi, ifreq)] = corr(t_score', tCoh(:,iroi,jroi,ifreq), 'Type','Spearman','tail','right');
            [rho_g(iroi,jroi,ifreq), pval_g(iroi,jroi, ifreq)] = corr(t_score', granger(:,iroi,jroi,ifreq), 'Type','Spearman');
            
        end
    end
end

[p_c_corrected, mask_c] = fdr(pval_c,0.05);
[p_g_corrected, mask_g] = fdr(pval_g,0.05);

%% sum across all freqs

figure
imagesc(sum(mask_c,3))
xlabel('Regions')
ylabel('Regions')
xticks=1:10;
xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
    'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',xticks,'YTickLabel',xticklabels)
xtickangle(45)
colorbar
title('MIM')

figure
imagesc(sum(mask_g,3))
xlabel('Regions')
ylabel('Regions')
xticks=1:10;
xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
    'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
set(gca,'XTick', xticks,'XTickLabel',xticklabels,'YTick',xticks,'YTickLabel',xticklabels)
xtickangle(45)
colorbar
title('GC') 

%% sum of significant freq bins per roi-roi-bin, divided by the number of freq bins in this band.
% So the colorbar is fixed from 0 to 1.

fbands = [1 4; 4 7; 7 15; 15 30; 30 45]; %theta, alpha, beta, low gamma, high gamma
fnames = {'theta', 'alpha', 'beta', 'low gamma', 'high gamma'};

figone(30,60)
for ibands = 1:size(fbands,1)
    subplot(2,5,ibands)
    imagesc(sum(mask_c(:,:,fbands(ibands,1):fbands(ibands,2)),3)./length(fbands(ibands,1):fbands(ibands,2)))
    title(['Corr Coh-TWSTRS ', fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45)
    caxis([0 1])
    
    subplot(2,5,ibands+5)
    imagesc(sum(mask_g(:,:,fbands(ibands,1):fbands(ibands,2)),3)./length(fbands(ibands,1):fbands(ibands,2)))
    title(['Corr GC-TWSTRS ' , fnames{ibands}])
    xlabel('Regions')
    ylabel('Regions')
    xticks=1:10;
    xticklabels = {'Precentral left','Precentral right', 'SMA left', 'SMA right', 'Parietal left',...
        'Parietal right', 'cerebellum','pallidum','LFP right','LFP left'};
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    xtickangle(45)
    caxis([0 1])
end

