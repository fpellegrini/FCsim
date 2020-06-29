function [mic, mim] = fp_get_mim_pca(A,CS,fqA,D,mode)
%mode is either a number for fixed pcs or 'max' (select npcs = rank of
%region data) or 'percent' (select npcs that 90% of the variance is
%preserved)

if isnumeric(mode)
    %fixed number of pcs for every roi 
    npcs = repmat(mode,D.nroi,1);
end

[nmeg, ni, nvox, nfreq] = size(A);

for aroi = 1:D.nroi
    
    %filter at current roi 
    clear A_ CSv
    A_ = A(:, :,D.ind_roi_cortex{aroi},:);
    nvoxroi = size(A_,3); %voxels in the current roi 
    A2{aroi} = reshape(A_, [nmeg, ni*nvoxroi, nfreq]);
    
    %project CS to voxel space
    for ifq = 1: nfreq
        CSv(:,:,ifq) = squeeze(A2{aroi}(:,:,fqA(ifq)))' * CS(:,:,ifq)...
            * squeeze(A2{aroi}(:,:,fqA(ifq)));
    end
    
    
    %zscoring
    clear CSz
    ZS{aroi} = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
    for ifreq = 1:nfreq
        CSz(:,:, ifreq) = ZS{aroi}'*squeeze(CSv(:, :,ifreq))*ZS{aroi};
    end
    
    
    %PCA
    clear CSs v v5 in V_ D_
    CSs = squeeze(sum(real(CSz),3)); %covariance
    [V_, D_] = eig(CSs);
    [D_, in] = sort(real(diag(D_)), 'descend');
    
    if strcmp(mode,'max')
        %pipeline 6)
        npcs(aroi) = min(rank(CSs),size(V_,1));
        
    elseif strcmp(mode,'percent')
        %pipeline 7)
        
        % variance explained
        vx_ = cumsum(D_)./sum(D_);
        % invx = 1:min(length(vx_), nmeg);
        npcs(aroi) = min(find(vx_>0.9));
    end
    
    V{aroi} = V_(:,in);
end

if strcmp(mode,'all')
    
    for ifi = 1:5
        npcs = repmat(ifi,D.nroi,1);
        [mic_fixed{ifi},mim_fixed{ifi}] = compute_mode_mim(ifi, D, npcs, V, A2, ZS, CS,fqA);
    end
    [mic_max,mim_max] = compute_mode_mim('max', D, npcs, V, A2, ZS, CS,fqA);
    [mic90,mim90] = compute_mode_mim('percent', D, npcs, V, A2, ZS, CS,fqA);
    
    mic.fixed = mic_fixed; 
    mic.max = mic_max; 
    mic.percent = mic90; 
    mim.fixed = mim_fixed; 
    mim.max = mim_max; 
    mim.percent = mim90; 
    
else
    [mic,mim] = compute_mode_mim(mode, D, npcs, V, A2, ZS, CS,fqA);
end

end


function [mic,mim] = compute_mode_mim(mode, D, npcs, V, A2, ZS, CS,fqA)

[nmeg, ~, nfreq] = size(A2{1});

%makes sure that rois have the same npcs in both hemispheres
if ~isnumeric(mode)
    
    partner_rois(1,:) = 1:D.nroi;
    partner_rois(2,[1:2:D.nroi-1])=2:2:D.nroi;
    partner_rois(2,[2:2:D.nroi]) = 1:2:D.nroi-1;
    
    if strcmp(mode,'max')
        for iroi = 1:D.nroi
            npcs(iroi) = min(npcs(iroi),npcs(partner_rois(2,iroi)));
        end
        
    elseif strcmp(mode,'percent')
        for iroi = 1:D.nroi
            npcs(iroi) = max(npcs(iroi),npcs(partner_rois(2,iroi)));
        end
    end
    
end


%concatenate filters
croi = 1;
for jroi = 1:D.nroi
    V{jroi} = V{jroi}(:, 1:npcs(jroi));
    for ifq = 1:nfreq
        P(:, croi:croi+npcs(jroi)-1,ifq) = A2{jroi}(:,:,ifq) * ZS{jroi} * real(V{jroi});
    end
    croi = croi +npcs(jroi);
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


%MIC and MIM
[mic,mim] =  fp_mim(Cohroi,npcs);

end


