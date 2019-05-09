%% group 

clear all
load('p_cluster_g_freq_abs.mat')

for iclus = 1:numel(p)
    if p(iclus) < 0.01
        clear a b c
        a = true_clu;
        a(a~=iclus) = 0;
        
        %spatial localization 
        b = sum(a,1);
        outname = sprintf('cluster_g_freq_abs_%d.nii',iclus);
        fp_data2nii(b,nan,[],outname)
        
        %spectral localization 
        c = sum(a,2);
        figure
        bar(c)
        xlabel('freqs')
        xticklabels = 0:5:92;
        xticks = linspace(1,length(c), numel(xticklabels));
        set(gca,'XTick', xticks,'XTickLabel',xticklabels)
        
        outname1 = sprintf('cluster_g_freq_abs_%d.png',iclus);
        print(outname1,'-dpng');
        close all

    end
   
end

%% singlesub

clear all
load('p_cluster_ss_c_freq_abs.mat')
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors


for isub = 1:numel(p)
    
    clear p1 true_clu
    p1 = p{isub};
    true_clu = TRUE_CLU{isub};
    
    for iclus = 1:numel(p1)
        if p1(iclus) < 0.05
            clear a b c
            a = true_clu;
            a(a~=iclus) = 0;

            %spatial localization 
            b = sum(a,1);
            outname = sprintf('cluster_ss_c_freq_abs_%d_Patient%s.nii',iclus,patientID{isub});
            fp_data2nii(b,nan,[],outname)

            %spectral localization 
            c = sum(a,2);
            figure
            bar(c)
            xlabel('freqs')
            xticklabels = 0:5:92;
            xticks = linspace(1,length(c), numel(xticklabels));
            set(gca,'XTick', xticks,'XTickLabel',xticklabels)

            outname1 = sprintf('cluster_ss_c_freq_abs_%dPatient%s.png',iclus,patientID{isub});
            print(outname1,'-dpng');
            close all

        end
    end
    
end

    
