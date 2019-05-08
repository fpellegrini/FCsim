load('p_cluster_g_freq_imag.mat')

for iclus = 1:numel(p)
    if p(iclus) < 0.01
        clear a b c
        a = true_clu;
        a(a~=iclus) = 0;
        
        %spatial localization 
        b = sum(a,1);
        outname = sprintf('a%d.nii',iclus);
        fp_data2nii(b,nan,[],outname)
        
        %spectral localization 
        c = sum(a,2);
        figure
        bar(c)
        xlabel('freqs')
        xticklabels = 0:5:92;
        xticks = linspace(1,length(c), numel(xticklabels));
        set(gca,'XTick', xticks,'XTickLabel',xticklabels)

        
    end
   
end


    
