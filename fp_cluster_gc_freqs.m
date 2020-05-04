function fp_cluster_gc(DIROUT, alpha,fwf,j)

%Group statistics on megmeg data.
%Clustering group statistics across space and frequencies.
%
%Input:
%
%abs_imag = either 'abs' or 'imag'
%
%testmethod = either 's' for signrank testing or 't' for simple
%thresholding
%
%fwf: testing method (less conservative): test first cluster against
%
%first, second with second etc. By default, clusters are always compared
%against the largest shuffled cluster.
%
%j=1 : only Julian's subjects


% fp_addpath_sabzi

if ~exist(DIROUT); mkdir(DIROUT); end

alpha_s = num2str(alpha);
alpha_s(1:2)=[];

if fwf==0
    fwf_s = [];
elseif fwf ==1
    fwf_s = 'fwf';
else
    error('Wrong fwf input')
end

if j == 0
    j_s = 'allsubs';
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
    [~, voxID] = fp_find_commonvox;
elseif j == 1
    j_s = 'j';
    patientID = {'04'; '07'; '09'; '10';'11';'20';'22';'25'};
    [pos, voxID] = fp_find_commonvox;
    voxID([3,7,8])=[];
else
    error('Wrong input j')
end

%load rois
nvox = size(pos,1);
for ii = 1: nvox
    [label{ii},code{ii},roi_id(ii)]=fp_get_mni_anatomy_new(pos(ii,:));
end
u_roi_id = sort(unique(roi_id));
nroi = numel(u_roi_id)-1;

%load GC
load(sprintf('%sDIFFGC',DIROUT));
[nsubs,nvox,nside,nfreq] = size(DIFFGC);
DIFFGC_s =DIFFGC;
nit = 1000;

%% freq wise all ROIs
DIFFGC= (sum(DIFFGC_s,2));
% true
fprintf('Testing...\n')
tic
[true_p,onoff,true_val,effectdir] = fp_get_signrank_results_gc(DIFFGC,alpha);
toc

% shuffled

for iit = 1:nit
    
    fprintf('Working on iteration %d \n',iit)
    clear onoff c_DIFFGC effectdir p_onoff n_onoff
    
    %shuffle signs of DIFFGC
    c_DIFFGC = DIFFGC .* (sign(randn(size(DIFFGC))));
    
    fprintf('Testing...\n')
    [~,~,shuf_val(iit,:,:,:), effectdir] = fp_get_signrank_results_gc(c_DIFFGC,alpha);
    
end
%
for iside=1:nside
    for ifq = 1:nfreq
        p(ifq,iside) = sum(squeeze(shuf_val(:,:,iside,ifq))> squeeze(true_val(:,iside,ifq)))/numel(shuf_val(:,:,iside,ifq));
    
    end

end 

for iside=1:2
    figure
    bar(squeeze(-log10(p(:,iside))))
    ylabel('-log10(p)')
    xlabel('Freqs')
    xticklabels = 0:5:92;
    xticks = linspace(1,size(p,1), numel(xticklabels));
    set(gca,'XTick', xticks,'XTickLabel',xticklabels)
    title(['side ' num2str(iside)])
end

%% freq-wise specific rois 

for aroi = 1:nroi
    
    clear true_val shuf_val p_roi
    DIFFGC= mean(sum(DIFFGC_s(:,roi_id==aroi,:,:),2));
    % true
    fprintf('Testing...\n')
    tic
    [~,~,true_val,~] = fp_get_signrank_results_gc(DIFFGC,alpha);
    toc
    
    % shuffled
    
    for iit = 1:nit
        clear onoff c_DIFFGC effectdir p_onoff n_onoff
        
        %shuffle signs of DIFFGC
        c_DIFFGC = DIFFGC .* (sign(randn(size(DIFFGC))));
        
        [~,~,shuf_val(iit,:,:,:), ~] = fp_get_signrank_results_gc(c_DIFFGC,alpha);
        
    end
    %
    for iside=1:nside
        for ifq = 1:nfreq
            p_roi(ifq,iside) = sum(squeeze(shuf_val(:,:,iside,ifq))> squeeze(true_val(:,iside,ifq)))/numel(shuf_val(:,:,iside,ifq));
            
        end
        
    end
    
    p_roi(p_roi==0)=0.00001;
    
    for iside=1:2
        figure
        bar(squeeze(-log10(p_roi(:,iside))))
        ylabel('-log10(p)')
        xlabel('Freqs')
        xticklabels = 0:5:92;
        xticks = linspace(1,size(p_roi,1), numel(xticklabels));
        set(gca,'XTick', xticks,'XTickLabel',xticklabels)
        title(['side ' num2str(iside), 'roi ' ROI_new.label{aroi}])
        
        outname = sprintf('%s_side%d.png', ROI_new.label{aroi},iside);
        print(outname,'-dpng');
        close all
    end
end

