% clear all
load('TWSTRS_gc.mat')
[pos, ~] = fp_find_commonvox;
sides = containers.Map([1 2], {'r', 'l'});


for iside = 1:2
    for iclus = 1:numel(p_neg{iside})
        if p_neg{iside}(iclus) < 0.05
            clear a b c
            %         a = true_clu;
            a = abs(true_rho(:,:,iside));
            a(true_clu_neg(:,:,iside)~=iclus) = 0;
            
            %spatial localization
            b = squeeze(mean(a(:,:),2));
%             figure; 
%             plot(b);
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

%%

for iside = 1:2
    for iclus = 1:numel(p_pos{iside})
        if p_pos{iside}(iclus) < 0.05
            clear a b c
            %         a = true_clu;
            a = abs(true_rho(:,:,iside));
            a(true_clu_pos(:,:,iside)~=iclus) = 0;
            
            %spatial localization
            b = squeeze(mean(a(:,:),2));
%             figure; 
%             plot(b);
            outname = sprintf('TWSTRS_gc_lcmv_%spos_%d.nii',sides(iside),iclus);
            fp_data2nii(b,pos,[],outname)
            
            %spectral localization
            c = sum(a,1);
            figure
            bar(c)
            xlabel('freqs')
            xticklabels = 0:5:92;
            xticks = linspace(1,length(c), numel(xticklabels));
            set(gca,'XTick', xticks,'XTickLabel',xticklabels)
            
            outname1 = sprintf('TWSTRS_gc_lcmv_%spos_%d.png',sides(iside), iclus);
            print(outname1,'-dpng');
            close all
            
        end
    end
    
end

%%

iside=2;
iclus = 1;

a = (true_rho(:,:,iside));
a(true_clu_neg(:,:,iside)~=iclus) = 0;

[iind,jind]= find(a==min(a(:)));

sub_gc = squeeze(diffgc_t(:,iind,1,jind));

corr(sub_gc,t_score')
scatter(t_score,sub_gc)
xlabel('twstrs score')
ylabel('GC (negative = cortex leads)')
outname = './lneg1_scatter_twstrs.png';
print(outname,'-dpng');
close all 


%% scatter 


iside=1;
iclus = 3;

a = (true_rho(:,:,iside));
a(true_clu_pos(:,:,iside)~=iclus) = 0;

[iind,jind]= find(a==max(a(:)));

sub_gc = squeeze(diffgc_t(:,iind,1,jind));

corr(sub_gc,t_score')
scatter(t_score,sub_gc)
xlabel('twstrs score')
ylabel('GC (negative = cortex leads)')

outname = './rpos3_scatter_twstrs.png';
print(outname,'-dpng');
close all 