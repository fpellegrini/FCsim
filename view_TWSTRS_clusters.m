clear all
load('TWSTRS_gc_lcmv.mat')
[pos, ~] = fp_find_commonvox;
sides = containers.Map([1 2], {'r', 'l'});


for iside = 1:2
    for iclus = 1:numel(p_neg{iside})
        if p_neg{iside}(iclus) < 0.01
            clear a b c
            %         a = true_clu;
            a = abs(true_rho(:,:,iside));
            a(true_clu_neg(:,:,iside)~=iclus) = 0;
            
            %spatial localization
            b = squeeze(mean(a(:,:),2));
            outname = sprintf('TWSTRS_gc_lcmv_%sneg_%d.nii',sides(iside),iclus);
            fp_data2nii(b,pos,[],outname)
            
            %spectral localization
            c = sum(a,1);
            figure
            bar(c)
            xlabel('freqs')
            xticklabels = 0:5:92;
            xticks = linspace(1,length(c), numel(xticklabels));
            set(gca,'XTick', xticks,'XTickLabel',xticklabels)
            
            outname1 = sprintf('TWSTRS_gc_lcmv_%sneg_%d.png',sides(iside), iclus);
            print(outname1,'-dpng');
            close all
            
        end
    end
    
end