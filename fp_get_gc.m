function [TRGC, GC,DIFFGC,to_save,t] = fp_get_gc(A,data, D,mode1,zs,t)
%mode1 is either a number for fixed pcs, or 'max' (select npcs = rank of
%region data), or 'percent' (select npcs that 90% of the variance is
%preserved), or 'case2' (only to pool dimensions, then summation), or
%'baseline', or 'all'
[n_sensors, ni, nvox] = size(A); 

tic
fprintf('Working on first part of gc_pca. \n')
%%
for aroi = 1:2%D.nroi
    
    %filter at current roi
    clear A_ datav
    A_ = A(:, :,D.ind_roi_cortex{aroi});
    nvoxroi(aroi) = size(A_,3); %voxels in the current roi
    A2{aroi} = reshape(A_, [n_sensors, ni*nvoxroi(aroi)]);
 
    %project to source space
    datav = A2{aroi}' * data(:,:);
        
   if ~strcmp(mode1,'case2')&& ~strcmp(mode1,'baseline')&& ~strcmp(mode1,'bandc')
        
        %zscoring
        clear CSz
        if zs
            dataz = zscore(datav(:,:)');
        else 
            dataz = datav';
        end
        
        %SVD
        clear data_ S_
        [data_, S_,~] = svds(double(dataz(:, :)),nvoxroi(aroi)*ni);
        
        % variance explained
        vx_ = cumsum(diag(S_).^2)./sum(diag(S_).^2);
        
        
        %npcs
        if isnumeric(mode1)
            %fixed number of pcs for every roi
            npcs(aroi) = mode1;
            var_explained(aroi) = vx_(npcs(aroi));
            
        elseif strcmp(mode1,'max')
            %pipeline 6)
            npcs(aroi) = min(find(vx_>0.99));
            
        elseif strcmp(mode1,'percent')
            %pipeline 7)
            npcs(aroi) = min(find(vx_>0.9));
            
%         elseif strcmp(mode1,'all')
%             
%             npcs.max(aroi) = min(find(vx_>0.99));           
%             npcs.percent(aroi) = min(find(vx_>0.9));
%             for ii = 1:5 
%                 var_explained(ii) = vx_(ii);
%             end
%             
        end
        
        % keep nPCA components with correct unit and scaling
       dataroi{aroi} = data_(:, 1:npcs(aroi))*S_(1:npcs(aroi), 1:npcs(aroi));
        
       
    else
        npcs=[];
        dataroi{aroi} = datav';
    end
end
t.pca = toc;
%%
fprintf('Working on compute_mode. \n')
if strcmp(mode1,'all')
    
    fprintf('fixed 1 to 5 \n')
    tic
    for ifi = 1:5
        npcs.fixed = repmat(ifi,D.nroi,1);
        [TRGC_fixed{ifi},GC_fixed{ifi},DIFFGC_fixed{ifi},to_save_fixed{ifi}] = ...
            fp_compute_mode_gc(ifi, D, npcs.fixed, V, A2, ZS, data_);
    end
    t.fixed = toc;
    
    fprintf('max \n')
    tic
    [TRGC_max,GC_max,DIFFGC_max,to_save_max] = fp_compute_mode_gc('max', D, npcs.max, V, A2, ZS, data_);
    t.ninetynine = toc;
    
    fprintf('90 percent \n')
    tic
    [TRGC90,GC90,DIFFGC90,to_save90] = fp_compute_mode_gc('percent', D, npcs.percent, V, A2, ZS, CS,data_);
    t.ninety = toc;
    
    fprintf('case2 and baseline \n')
    tic
    [TRGC_bandc,GC_bandc,DIFFGC_bandc, to_save_bandc] = fp_compute_mode_gc('bandc', D, [], [], A2, [], datav);
    t.bandc = toc;
    
    TRGC.fixed = TRGC_fixed;
    TRGC.max = TRGC_max;
    TRGC.percent = TRGC90;
    TRGC.case2 = TRGC_bandc.case2;
    TRGC.baseline = TRGC_bandc.baseline;
    GC.fixed = GC_fixed;
    GC.max = GC_max;
    GC.percent = GC90;
    GC.case2 = GC_bandc.case2;
    GC.baseline = GC_bandc.baseline;
    DIFFGC.fixed = DIFFGC_fixed;
    DIFFGC.max = DIFFGC_max;
    DIFFGC.percent = DIFFGC90;
    DIFFGC.case2 = DIFFGC_bandc.case2;
    DIFFGC.baseline = DIFFGC_bandc.baseline;

    to_save.fixed = to_save_fixed;
    for ii = 1:5
        to_save.fixed{ii}.var_explained = var_explained(ii);
    end
    to_save.max = to_save_max;
    to_save.percent = to_save90;
    to_save.bandc = to_save_bandc; %too large to be saved?
    to_save.nvoxroi = nvoxroi;
    
    %% Correlations 
   tic
    nvoxroi_all = nvoxroi'*nvoxroi;
    nvoxroi_all = nvoxroi_all(:);
    for ii = 1:5
        c1 = sum(TRGC.fixed{ii},3);
        c2 = sum(GC.fixed{ii},3);
        c3 = sum(DIFFGC.fixed{ii},3);
        to_save.fixed{ii}.corr_voxtrgc = corr(nvoxroi_all,c1(:));
        to_save.fixed{ii}.corr_voxgc = corr(nvoxroi_all,c2(:));
        to_save.fixed{ii}.corr_voxdiffgc = corr(nvoxroi_all,c3(:));    
        
    end
    c1 = sum(TRGC.max,3);
    c2 = sum(GC.max,3);
    c3 = sum(DIFFGC.max,3);
    to_save.max.corr_voxtrgc = corr(nvoxroi_all,c1(:));
    to_save.max.corr_voxgc = corr(nvoxroi_all ,c2(:));
    to_save.max.corr_voxnpcs = corr(nvoxroi', to_save.max.npcs');
    to_save.max.corr_voxdiffgc = corr(nvoxroi_all,c3(:));
                
    c1 = sum(TRGC.percent,3);
    c2 = sum(GC.percent,3);
    c3 = sum(DIFFGC.percent,3);
    to_save.percent.corr_voxtrgc = corr(nvoxroi_all,c1(:));
    to_save.percent.corr_voxgc = corr(nvoxroi_all,c2(:));
    to_save.percent.corr_voxnpcs = corr(nvoxroi', to_save.percent.npcs');
    to_save.percent.corr_voxdiffgc = corr(nvoxroi_all,c3(:));
            
    c1 = sum(TRGC.case2,3);
    c2 = sum(GC.case2,3);
    c3 = sum(DIFFGC.case2,3);
    to_save.case2.corr_voxtrgc = corr(nvoxroi_all,c1(:));
    to_save.case2.corr_voxgc = corr(nvoxroi_all,c2(:));
    to_save.case3.corr_voxdiffgc = corr(nvoxroi_all,c3(:));
    
    c1 = sum(TRGC.baseline,3);
    c2 = sum(GC.baseline,3);
    c3 = sum(DIFFGC.baseline,3);
    to_save.baseline.corr_voxtrgc = corr(nvoxroi_all,c1(:));
    to_save.baseline.corr_voxgc = corr(nvoxroi_all,c2(:));
    to_save.baseline.corr_voxdiffgc = corr(nvoxroi_all,c3(:));

    t.corrs = toc;
    
    
else
    tic
    [TRGC,GC,DIFFGC, to_save] = fp_compute_mode_gc(mode1, D, npcs, V, A2, ZS, CS,fqA,nfqA);
    t.gc = toc;
end






