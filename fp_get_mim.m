function [mic, mim,to_save] = fp_get_mim(A,CS,fqA,nfqA, D,ihemi,mode1)
%mode1 is either a number for fixed pcs, or 'max' (select npcs = rank of
%region data), or 'percent' (select npcs that 90% of the variance is
%preserved), or 'case2' (mim only to pool dimensions, then summation), or
%'baseline', or 'all'
[nmeg, ni, nvox,~] = size(A);
nfreq = size(CS,3); 

fprintf('Working on first part of mim_pca. \n')
tic
for aroi = 1:D.nroi
    
    %filter at current roi
    clear A_ CSv
    A_ = A(:, :,D.ind_roi_cortex{aroi},:);
    nvoxroi = size(A_,3); %voxels in the current roi
    A2{aroi} = reshape(A_, [nmeg, ni*nvoxroi, nfqA]);
    
    
    if ~strcmp(mode1,'case2')&& ~strcmp(mode1,'baseline')
        
        %project CS to voxel space
        for ifq = 1: nfreq
            CSv(:,:,ifq) = squeeze(A2{aroi}(:,:,fqA(ifq)))' * CS(:,:,ifq)...
                * squeeze(A2{aroi}(:,:,fqA(ifq)));
        end
        
        %zscoring
        clear CSz
        ZS{aroi} = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
        
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
        
        %npcs
        if isnumeric(mode1)
            %fixed number of pcs for every roi
            npcs(aroi) = mode1;
            
        elseif strcmp(mode1,'max')
            %pipeline 6)
            
            % npcs(aroi) = min(nmeg,rank(CSs));
            vx_ = cumsum(D_)./sum(D_);
            npcs(aroi) = min(find(vx_>0.99999));
            
        elseif strcmp(mode1,'percent')
            %pipeline 7)
            
            % variance explained
            vx_ = cumsum(D_)./sum(D_);
            npcs(aroi) = min(find(vx_>0.9));
            
        elseif strcmp(mode1,'all')
            npcs.max(aroi) = min(rank(CSs),size(V_,1));
            vx_ = cumsum(D_)./sum(D_);
            npcs.percent(aroi) = min(find(vx_>0.9));
            
        end
    else
        
        ZS = [];
        V=[];
        npcs=[];
    end
end
toc
%%
fprintf('Working on compute_mode. \n')
if strcmp(mode1,'all')
    
    fprintf('fixed 1 to 5 \n')
    tic
    for ifi = 1:4
        npcs.fixed = repmat(ifi,D.nroi,1);
        [mic_fixed{ifi},mim_fixed{ifi},to_save_fixed{ifi}] = fp_compute_mode_mim(ifi, D, npcs.fixed, V, A2, ZS, CS,fqA,nfqA,ihemi);
    end
    toc
    
    fprintf('max \n')
    tic
    [mic_max,mim_max,to_save_max] = fp_compute_mode_mim('max', D, npcs.max, V, A2, ZS, CS,fqA,nfqA,ihemi);
    toc
    
    fprintf('90 percent \n')
    tic
    [mic90,mim90,to_save90] = fp_compute_mode_mim('percent', D, npcs.percent, V, A2, ZS, CS,fqA,nfqA,ihemi);
    toc
    
    fprintf('case2 and baseline \n')
    tic
    [mic_bandc,mim_bandc,to_save_bandc] = fp_compute_mode_mim('bandc',D,[],[],A2,[],CS,fqA,nfqA,ihemi);
    toc
    
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
    to_save.max = to_save_max;
    to_save.percent = to_save90;
%     to_save.bandc = to_save_bandc; %too large to be saved 
    
    
else
    [mic,mim,to_save] = fp_compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA,nfqA,ihemi);
end





