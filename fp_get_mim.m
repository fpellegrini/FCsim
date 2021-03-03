function [mic, mim,to_save, mean_icoh, mean_acoh,t] = fp_get_mim(A,CS,fqA,nfqA, D,ihemi,mode1,zs,t)
%mode1 is either a number for fixed pcs, or 'max' (select npcs = rank of
%region data), or 'percent' (select npcs that 90% of the variance is
%preserved), or 'case2' (mim only to pool dimensions, then summation), or
%'baseline', or 'all'
[nmeg, ni, nvox,~] = size(A);
nfreq = size(CS,3); 

tic
fprintf('Working on first part of mim_pca. \n')
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
            for ii = 1:6 
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
    
    fprintf('fixed 1 to 6 \n')
    tic
    for ifi = 1:6
        npcs.fixed = repmat(ifi,D.nroi,1);
        [mic_fixed{ifi},mim_fixed{ifi},to_save_fixed{ifi},mean_icoh_fixed{ifi}, mean_acoh_fixed{ifi}] = fp_compute_mode_mim(ifi, D, npcs.fixed, V, A2, ZS, CS,fqA,nfqA,ihemi);
    end
    t.fixed = toc;
    
    fprintf('max \n')
    tic
    [mic_max,mim_max,to_save_max,mean_icoh_max,mean_acoh_max] = fp_compute_mode_mim('max', D, npcs.max, V, A2, ZS, CS,fqA,nfqA,ihemi);
    t.ninetynine = toc;
    
    fprintf('90 percent \n')
    tic
    [mic90,mim90,to_save90,mean_icoh_90, mean_acoh_90] = fp_compute_mode_mim('percent', D, npcs.percent, V, A2, ZS, CS,fqA,nfqA,ihemi);
    t.ninety = toc;
    
    fprintf('case2 and baseline \n')
    tic
    [mic_bandc,mim_bandc,to_save_bandc,~] = fp_compute_mode_mim('bandc',D,[],[],A2,[],CS,fqA,nfqA,ihemi);
    t.bandc = toc;
    
    mic.fixed = mic_fixed;
    mic.max = mic_max;
    mic.percent = mic90;
    mic.case2 = mic_bandc.case2;
    mic.baseline = mic_bandc.baseline;
    mim.fixed = mim_fixed;
    mim.max = mim_max;
    mim.percent = mim90;
    mim.case2 = mim_bandc.case2;
    mim.baseline = mim_bandc.baseline;

    to_save.fixed = to_save_fixed;
    for ii = 1:6
        to_save.fixed{ii}.var_explained = var_explained(ii);
    end
    to_save.max = to_save_max;
    to_save.percent = to_save90;
    to_save.bandc = to_save_bandc; %too large to be saved?
    to_save.nvoxroi = nvoxroi;
    
    
    mean_icoh.fixed = mean_icoh_fixed; 
    mean_icoh.max = mean_icoh_max; 
    mean_icoh.percent = mean_icoh_90; 
    
    mean_acoh.fixed = mean_acoh_fixed; 
    mean_acoh.max = mean_acoh_max; 
    mean_acoh.percent = mean_acoh_90; 
    
    %% Correlations 
    tic
    nvoxroi_all = nvoxroi'*nvoxroi;
    nvoxroi_all = nvoxroi_all(:);
    for ii = 1:6
        c1 = sum(mim.fixed{ii},3);
        c2 = sum(mic.fixed{ii},3);
        c3 = sum(mean_icoh.fixed{ii},3);
        c4 = sum(mean_acoh.fixed{ii},3);
        to_save.fixed{ii}.corr_voxmim = corr(nvoxroi_all,c1(:));
        to_save.fixed{ii}.corr_voxmic = corr(nvoxroi_all,c2(:));
        to_save.fixed{ii}.corr_voxmeancoh = corr(nvoxroi_all,c3(:));        
        to_save.fixed{ii}.corr_voxmeanabscoh = corr(nvoxroi_all,c4(:));
        
    end
    c1 = sum(mim.max,3);
    c2 = sum(mic.max,3);
    c3 = sum(mean_icoh.max,3);
    c4 = sum(mean_acoh.max,3);
    to_save.max.corr_voxmim = corr(nvoxroi_all,c1(:));
    to_save.max.corr_voxmic = corr(nvoxroi_all ,c2(:));
    to_save.max.corr_voxnpcs = corr(nvoxroi', to_save.max.npcs');
    to_save.max.corr_voxmeancoh = corr(nvoxroi_all,c3(:));
    to_save.max.corr_voxmeanabscoh = corr(nvoxroi_all,c4(:));
                
    c1 = sum(mim.percent,3);
    c2 = sum(mic.percent,3);
    c3 = sum(mean_icoh.percent,3);
    c4 = sum(mean_acoh.percent,3);
    to_save.percent.corr_voxmim = corr(nvoxroi_all,c1(:));
    to_save.percent.corr_voxmic = corr(nvoxroi_all,c2(:));
    to_save.percent.corr_voxnpcs = corr(nvoxroi', to_save.percent.npcs');
    to_save.percent.corr_voxmeancoh = corr(nvoxroi_all,c3(:));
    to_save.percent.corr_voxmeanabscoh = corr(nvoxroi_all,c4(:));
            
    c1 = sum(mim.case2,3);
    c2 = sum(mic.case2,3);
    to_save.case2.corr_voxmim = corr(nvoxroi_all,c1(:));
    to_save.case2.corr_voxmic = corr(nvoxroi_all,c2(:));
    
    c1 = sum(mim.baseline,3);
    c2 = sum(mic.baseline,3);
    to_save.baseline.corr_voxmim = corr(nvoxroi_all,c1(:));
    to_save.baseline.corr_voxmic = corr(nvoxroi_all,c2(:));

    t.corrs = toc;
    
    
else
    tic
    [mic,mim,to_save, mean_icoh, mean_acoh] = fp_compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA,nfqA,ihemi);
    t.mim = toc;
end







