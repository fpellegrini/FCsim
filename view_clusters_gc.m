clear all

load('p_gc_lcmv_01__allsubs.mat')
[pos, ~] = fp_find_commonvox;
sides = containers.Map([1 2], {'r', 'l'});
vals = containers.Map([1 2], {'pos', 'neg'});


for iside = 1:2 
    for ival = 1:2 
        for iclus = 1:numel(p{iside,ival})
            if p{iside,ival}(iclus) < 0.05
                clear a b c
                a = squeeze(true_clu(:,:,iside,ival))';
        %         a = true_testval';
                a(a~=iclus) = 0;

                %spatial localization
                b = squeeze(sum(a,1));
                outname = sprintf('cluster_gc_01_%s%s_%d_lcmv.nii',sides(iside), vals(ival),iclus);
                fp_data2nii(b,pos,[],outname)

                %spectral localization
                c = sum(a,2);
                figure
                bar(c)
                xlabel('freqs')
                xticklabels = 0:5:92;
                xticks = linspace(1,length(c), numel(xticklabels));
                set(gca,'XTick', xticks,'XTickLabel',xticklabels)

                outname1 = sprintf('cluster_gc_01_%s%s_%d_lcmv.png',sides(iside), vals(ival), iclus);
                print(outname1,'-dpng');
                close all

            end

        end
    end 
end 

%% clear all

load('TWSTRS_gc.mat')
[pos, ~] = fp_find_commonvox;
sides = containers.Map([1 2], {'r', 'l'});
vals = containers.Map([1 2], {'pos', 'neg'});


for iside = 1:2 
        for iclus = 1:numel(p_neg{iside})
            if p_neg{iside}(iclus) < 0.05
                clear a b c
                a = true_rho(:,:,iside);
                a(true_clu_neg(:,:,iside)~=iclus) = 0;

                %spatial localization
                b = squeeze(sum(a,2));
                outname = sprintf('TWSTRS_gc_neg_%s_%d.nii',sides(iside),iclus);
                fp_data2nii(b,pos,[],outname)

                %spectral localization
                c = sum(a,1);
                figure
                bar(c)
                xlabel('freqs')
                xticklabels = 0:5:92;
                xticks = linspace(1,length(c), numel(xticklabels));
                set(gca,'XTick', xticks,'XTickLabel',xticklabels)

                outname1 = sprintf('TWSTRS_gc_neg_%s_%d.png',sides(iside), iclus);
                print(outname1,'-dpng');
                close all

            end

        end 
end 


for iside = 1:2 
        for iclus = 1:numel(p_pos{iside})
            if p_pos{iside}(iclus) < 0.05
                clear a b c
                a = true_rho(:,:,iside);
                a(true_clu_pos(:,:,iside)~=iclus) = 0;

                %spatial localization
                b = squeeze(sum(a,2));
                outname = sprintf('TWSTRS_gc_pos_%s_%d.nii',sides(iside),iclus);
                fp_data2nii(b,pos,[],outname)

                %spectral localization
                c = sum(a,1);
                figure
                bar(c)
                xlabel('freqs')
                xticklabels = 0:5:92;
                xticks = linspace(1,length(c), numel(xticklabels));
                set(gca,'XTick', xticks,'XTickLabel',xticklabels)

                outname1 = sprintf('TWSTRS_gc_pos_%s_%d.png',sides(iside), iclus);
                print(outname1,'-dpng');
                close all

            end

        end 
end 