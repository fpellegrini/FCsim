%% group

clear all
load('p_cluster_g_c_freq_imag_test2D.mat')
[pos, ~] = fp_find_commonvox;

for iclus = 1:numel(p)
    if p(iclus) < 0.01
        clear a b c
        a = true_clu;
%         a = true_testval';
        a(true_clu~=iclus) = 0;
        
        %spatial localization
        b = squeeze(sum(a,1));
        outname = sprintf('cluster_g_c_freq_imag_test2D_%d.nii',iclus);
        fp_data2nii(b,pos,[],outname)
        
        %spectral localization
        c = sum(a,2);
        figure
        bar(c)
        xlabel('freqs')
        xticklabels = 0:5:92;
        xticks = linspace(1,length(c), numel(xticklabels));
        set(gca,'XTick', xticks,'XTickLabel',xticklabels)
        
        outname1 = sprintf('cluster_g_freq_imag_test2D_%d.png',iclus);
        print(outname1,'-dpng');
        close all
        
    end
    
end

%% group bands 
afs = {'theta','alpha','beta','gamma_low','gamma_high'};
abs_imag = {'abs','imag'};

for iafs = [1:numel(afs)]
    for iabs = 2
        clearvars -except iafs afs iabs abs_imag
        load(sprintf('p_cluster_g_%s_%s_j2.mat',afs{iafs},abs_imag{iabs}))

        for iclus = 1:numel(p)
            if p(iclus) < 0.01
                clear a b c
                a = true_clu;
                a(a~=iclus) = 0;

                %spatial localization
                b = sum(a,1);
                outname = sprintf('cluster_g_%s_%s_j2_%d.nii',afs{iafs},abs_imag{iabs},iclus);
                fp_data2nii(b,nan,[],outname)

            end

        end
    end 
end


%% singlesub

clear all
load('p_cluster_ss_c_freq_imag.mat')
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors


for isub = 1:numel(p)
    
    clear p1 true_clu
    p1 = p{isub};
    true_clu = TRUE_CLU{isub};
    
    if ~isempty(true_clu)
        for iclus = 1:numel(p1)
            if p1(iclus) < 0.05
                clear a b c
                a = true_clu;
                a(a~=iclus) = 0;
                
                %spatial localization
                b = sum(a,1);
                
                outname = sprintf('cluster_ss_c_freq_imag_%d_Patient%s.nii',iclus,patientID{isub});
                fp_data2nii(b,nan,[],outname)
                
                %spectral localization
                c = sum(a,2);
                figure
                bar(c)
                xlabel('freqs')
                xticklabels = 0:5:92;
                xticks = linspace(1,length(c), numel(xticklabels));
                set(gca,'XTick', xticks,'XTickLabel',xticklabels)
                
                outname1 = sprintf('cluster_ss_c_freq_imag_%dPatient%s.png',iclus,patientID{isub});
                print(outname1,'-dpng');
                close all
                
            end
        end
    end
    
    
end


