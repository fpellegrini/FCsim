function [mic, mim] = fp_get_mim(A,CS,fqA,D,mode1)
%mode1 is either a number for fixed pcs, or 'max' (select npcs = rank of
%region data), or 'percent' (select npcs that 90% of the variance is
%preserved), or 'case2' (mim only to pool dimensions, then summation), or
%'baseline', or 'all'

[nmeg, ni, nvox, nfreq] = size(A);

fprintf('Working on first part of mim_pca. \n')
tic
for aroi = 1:D.nroi
    
    %filter at current roi
    clear A_ CSv
    A_ = A(:, :,D.ind_roi_cortex{aroi},:);
    nvoxroi = size(A_,3); %voxels in the current roi
    A2{aroi} = reshape(A_, [nmeg, ni*nvoxroi, nfreq]);
    
    
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
            npcs(aroi) = min(rank(CSs),size(V_,1));
            
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

fprintf('Working on compute_mode. \n')
tic
if strcmp(mode1,'all')
    
    fprintf('fixed 1 to 5 \n')
    tic
    for ifi = 1:5
        npcs.fixed = repmat(ifi,D.nroi,1);
        [mic_fixed{ifi},mim_fixed{ifi}] = compute_mode_mim(ifi, D, npcs.fixed, V, A2, ZS, CS,fqA);
    end
    toc
    fprintf('max \n')
    tic
    [mic_max,mim_max] = compute_mode_mim('max', D, npcs.max, V, A2, ZS, CS,fqA);
    toc
    tic
    fprintf('90 percent \n')
    toc
    tic
    [mic90,mim90] = compute_mode_mim('percent', D, npcs.percent, V, A2, ZS, CS,fqA);
    toc
    tic
    fprintf('case2 and baseline \n')
    toc
    tic
    [mic_bandc,mim_bandc] = compute_mode_mim('bandc',D,[],[],A2,[],CS,fqA);
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
    
    
else
    [mic,mim] = compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA);
end

toc
end



function [mic,mim] = compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA)

[nmeg, ~, nfreq] = size(A2{1});

%makes sure that rois have the same npcs in both hemispheres
if strcmp(mode1,'max')||strcmp(mode1,'percent')
    
    partner_rois(1,:) = 1:D.nroi;
    partner_rois(2,[1:2:D.nroi-1])=2:2:D.nroi;
    partner_rois(2,[2:2:D.nroi]) = 1:2:D.nroi-1;
    
    if strcmp(mode1,'max')
        for iroi = 1:D.nroi
            npcs(iroi) = min(npcs(iroi),npcs(partner_rois(2,iroi)));
        end
        
    elseif strcmp(mode1,'percent')
        for iroi = 1:D.nroi
            npcs(iroi) = max(npcs(iroi),npcs(partner_rois(2,iroi)));
        end
    end
    
end

%concatenate filters
if strcmp(mode1,'bandc')|| strcmp(mode1,'baseline')||strcmp(mode1,'case2')
    P=[];
    for jroi = 1:D.nroi
        P = cat(2,P,A2{jroi});
    end
else
    croi = 1;
    for jroi = 1:D.nroi
        V{jroi} = V{jroi}(:, 1:npcs(jroi));
        for ifq = 1:nfreq
            P(:, croi:croi+npcs(jroi)-1,ifq) = A2{jroi}(:,:,ifq) * ZS{jroi} * real(V{jroi});
        end
        croi = croi +npcs(jroi);
    end
end


%apply all filters
CSroi = [];
for ifreq = 1:nfreq
    CSroi(:, :, ifreq) = reshape(P(:,:,fqA(ifreq)), nmeg, [])'*CS(:, :, ifreq)...
        *reshape(P(:,:,fqA(ifreq)), nmeg, []);
end


%divide by power to obtain coherence
clear Cohroi
for ifreq = 1: nfreq
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end
clear CSroi


if strcmp(mode1,'baseline')||strcmp(mode1,'bandc')
    
    baseline = Cohroi(D.sub_ind_cortex,D.sub_ind_cortex,:);
    iReg = length(D.sub_ind_roi{1});
    
    if iReg ==1 %in this case, mim and mic would not make sense 
        mic = baseline;
        mim=mic;
    else
        %mim across iReg voxels of one region
        [mic,mim]= fp_mim(baseline,npcs);
    end
    
    if strcmp(mode1,'bandc')
        a = mic; 
        b = mim; 
        clear mim mic
    end
end
    
if strcmp(mode1,'case2')|| strcmp(mode1,'bandc')
    
    %MIC and MIM
    [mic_v,mim_v] =  fp_mim(Cohroi,npcs);
    
    %sum up mim/ mic within rois    
    for iroi = 1:D.nroi
        for jroi = 1:D.nroi
            mic(iroi,jroi,:) = squeeze(sum(sum(mic_v(D.ind_roi_cortex{iroi},...
                D.ind_roi_cortex{jroi},:),1),2));
            mim(iroi,jroi,:) = squeeze(sum(sum(mim_v(D.ind_roi_cortex{iroi},...
                D.ind_roi_cortex{jroi},:),1),2));
        end
    end
    
    if strcmp(mode1,'bandc')
        d = mic; 
        e = mim; 
        clear mim mic
        mic.baseline = a;
        mim.baseline = b; 
        mic.case2 = d; 
        mim.case2 = e; 
    end
    
else
    
    %MIC and MIM
    [mic,mim] =  fp_mim(Cohroi,npcs);
end

end


