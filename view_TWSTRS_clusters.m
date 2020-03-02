clear all
load('TWSTRS.mat')
[pos, ~] = fp_find_commonvox;

for iclus = 1:numel(p_pos)
    if p_pos(iclus) < 0.01
        clear a b c
%         a = true_clu;
        a = true_rho;
        a(true_clu_pos~=iclus) = 0;
        
        %spatial localization
        b = squeeze(mean(a(4:13,:),1));
        outname = sprintf('TWSTRS_%d_alpha.nii',iclus);
        fp_data2nii(b,pos,[],outname)
        
        %spectral localization
        c = sum(a,2);
        figure
        bar(c)
        xlabel('freqs')
        xticklabels = 0:5:92;
        xticks = linspace(1,length(c), numel(xticklabels));
        set(gca,'XTick', xticks,'XTickLabel',xticklabels)
        
        outname1 = sprintf('TFCE_d_comp_%d.png',iclus);
        print(outname1,'-dpng');
        close all
        
    end
    
end