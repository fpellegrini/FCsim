clear all

load('p_gc_01__allsubs.mat')
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
                outname = sprintf('cluster_gc_01_%s%s_%d.nii',sides(iside), vals(ival),iclus);
                fp_data2nii(b,pos,[],outname)

                %spectral localization
                c = sum(a,2);
                figure
                bar(c)
                xlabel('freqs')
                xticklabels = 0:5:92;
                xticks = linspace(1,length(c), numel(xticklabels));
                set(gca,'XTick', xticks,'XTickLabel',xticklabels)

                outname1 = sprintf('cluster_gc_01_%s%s_%d.png',iclus);
                print(outname1,'-dpng');
                close all

            end

        end
    end 
end 