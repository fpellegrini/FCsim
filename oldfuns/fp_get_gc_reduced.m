function [TRGC,GC,DIFFGC,to_save,t] = fp_get_gc_reduced(A,CS,fqA,nfqA, D,mode1,zs,t)
%mode1 is either a number for fixed pcs, or 'max' (select npcs = rank of
%region data), or 'percent' (select npcs that 90% of the variance is
%preserved), or 'case2' (mvgc only to pool dimensions, then summation), or
%'baseline', or 'all'
[nmeg, ni, nvox,~] = size(A);
nfreq = size(CS,3); 

tic
fprintf('Working on first part of gc_pca. \n')
%%
for aroi = 1:D.nroi
    
    %filter at current roi
    clear A_ CSv
    A_ = A(:, :,D.ind_roi_cortex{aroi},:);
    nvoxroi(aroi) = size(A_,3); %voxels in the current roi
    A2{aroi} = reshape(A_, [nmeg, ni*nvoxroi(aroi), nfqA]);
    
    
    if ~strcmp(mode1,'case2')&& ~strcmp(mode1,'baseline')&& ~strcmp(mode1,'bandc')
        
        %project CS to voxel space
        for ifq = 1: nfreq
            CSv(:,:,ifq) = squeeze(A2{aroi}(:,:,fqA(ifq)))' * CS(:,:,ifq)...
                * squeeze(A2{aroi}(:,:,fqA(ifq)));
        end
        
        %zscoring
        clear CSz
        if zs
            ZS{aroi} = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
        else 
            ZS{aroi} = eye(size(CSv,1));
        end
        
        %apply ZS
        for ifreq = 1:nfreq
            CSz(:,:, ifreq) = ZS{aroi}'*squeeze(CSv(:, :,ifreq))*ZS{aroi};
        end
        
        %PCA
        clear CSs v v5 in V_ D_
        CSs = squeeze(sum(real(CSz),3)); %covariance
        [V_, D_] = eig(CSs);
        [D_, in] = sort(real(diag(D_)), 'descend');
        V{aroi} = V_(:,in);
        vx_ = cumsum(D_)./sum(D_);
        
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
            
        elseif strcmp(mode1,'all')
            
            npcs.max(aroi) = min(find(vx_>0.99));           
            npcs.percent(aroi) = min(find(vx_>0.9));
            for ii = 1:5 
                var_explained(ii) = vx_(ii);
            end
            
        end
    else
        
        ZS = [];
        V=[];
        npcs=[];
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
            fp_compute_mode_gc(ifi, D, npcs.fixed, V, A2, ZS, CS,fqA,nfqA);
    end
    t.fixed = toc;
    
    fprintf('max \n')
    tic
    [TRGC_max,GC_max,DIFFGC_max,to_save_max] = fp_compute_mode_gc('max', D, npcs.max, V, A2, ZS, CS,fqA,nfqA);
    t.ninetynine = toc;
    
    fprintf('90 percent \n')
    tic
    [TRGC90,GC90,DIFFGC90,to_save90] = fp_compute_mode_gc('percent', D, npcs.percent, V, A2, ZS, CS,fqA,nfqA);
    t.ninety = toc;
    
    fprintf('baseline \n')
    tic
    [TRGC_b,GC_b,DIFFGC_b, to_save_b] = fp_compute_mode_gc('baseline', D, [], [], A2, [], CS,fqA,nfqA);
    t.bandc = toc;
    
    TRGC.fixed = TRGC_fixed;
    TRGC.max = TRGC_max;
    TRGC.percent = TRGC90;
    TRGC.baseline = TRGC_b;
    GC.fixed = GC_fixed;
    GC.max = GC_max;
    GC.percent = GC90;
    GC.baseline = GC_b;
    DIFFGC.fixed = DIFFGC_fixed;
    DIFFGC.max = DIFFGC_max;
    DIFFGC.percent = DIFFGC90;
    DIFFGC.baseline = DIFFGC_b;
    
    to_save.fixed = to_save_fixed;
    for ii = 1:5
        to_save.fixed{ii}.var_explained = var_explained(ii);
    end
    to_save.max = to_save_max;
    to_save.percent = to_save90;
    to_save.bandc = to_save_b; %too large to be saved?
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







