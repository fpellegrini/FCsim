function [mic, mim] = fp_get_mim(A,CS,fqA,D,ihemi,mode1)
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
            npcs(aroi) = min(nmeg,size(V_,1));
            
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
if strcmp(mode1,'all')
    
    fprintf('fixed 1 to 5 \n')
    tic
    for ifi = 1:5
        npcs.fixed = repmat(ifi,D.nroi,1);
        [mic_fixed{ifi},mim_fixed{ifi}] = compute_mode_mim(ifi, D, npcs.fixed, V, A2, ZS, CS,fqA,ihemi);
    end
    toc
    
    fprintf('max \n')
    tic
    [mic_max,mim_max] = compute_mode_mim('max', D, npcs.max, V, A2, ZS, CS,fqA,ihemi);
    toc
    
    fprintf('90 percent \n')
    tic
    [mic90,mim90] = compute_mode_mim('percent', D, npcs.percent, V, A2, ZS, CS,fqA,ihemi);
    toc
    
    fprintf('case2 and baseline \n')
    tic
    [mic_bandc,mim_bandc] = compute_mode_mim('bandc',D,[],[],A2,[],CS,fqA,ihemi);
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
    [mic,mim] = compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA,ihemi);
end

end



function [mic,mim] = compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA,ihemi)

ndim = 2; 

%when ihemi == 1: makes sure that rois have the same npcs in both hemispheres
if (strcmp(mode1,'max')||strcmp(mode1,'percent')) && (ihemi==1)    
    npcs = fp_symmetrize_hemispheres(D,npcs,mode1);    
end

%loop over all roi combinations 
for oroi = 1:D.nroi
    for uroi = oroi+1:D.nroi
        
        clear P
        [nmeg, dummy, nfreq] = size(A2{oroi});
        nvoxreg1 = dummy/ndim;
     
        %concatenate filters
        if strcmp(mode1,'bandc')|| strcmp(mode1,'baseline')||strcmp(mode1,'case2')
            P = cat(2,A2{oroi},A2{uroi});
            
        else
            croi = 1;
            for jroi = [oroi uroi]
                for ifq = 1:nfreq
                    P(:, croi:croi+npcs(jroi)-1,ifq) = A2{jroi}(:,:,ifq) * ZS{jroi} * real(V{jroi}(:, 1:npcs(jroi)));
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
            %preselect "true" voxels (ind sub_ind_roi_region is the
            %index within the roi. Caution with 2 dimensions! 
            
            %indeces of the true voxels in Cohroi
            clear ouind1 ouind2 baseline npcs 
            ouind1= cat(2,(D.sub_ind_roi_region{oroi}.*ndim)-1,(D.sub_ind_roi_region{oroi}.*ndim));
            ouind2= cat(2,(D.sub_ind_roi_region{uroi}.*ndim)-1+(nvoxreg1*ndim),...
                (D.sub_ind_roi_region{uroi}.*ndim)+(nvoxreg1*ndim));
            
            baseline = Cohroi([ouind1 ouind2],[ouind1 ouind2],:);    
            iReg = length(D.sub_ind_roi{1});

            %mim across iReg voxels and ndim of one region
            npcs = repmat(iReg*ndim,size(baseline,1)/(ndim*iReg),1);
            clear a b
            [a, b]= fp_mim(baseline,npcs);
            mic1([oroi uroi],[oroi uroi],:) = a;
            mim1([oroi uroi],[oroi uroi],:) = b;

            
        end
        
        if strcmp(mode1,'case2')|| strcmp(mode1,'bandc')
            
            %MIC and MIM
            npcs = repmat(ndim,size(Cohroi,1)/(ndim),1);
            [mic_v,mim_v] =  fp_mim(Cohroi,npcs);
            
            %sum up mim/ mic within rois
            mic2(oroi,uroi,:) = squeeze(mean(mean(mic_v(1:nvoxreg1,nvoxreg1+1:end,:),1),2));
            mic2(uroi,oroi,:) = mic2(oroi,uroi,:);
            mic2(oroi,oroi,:) = squeeze(mean(mean(mic_v(1:nvoxreg1,1:nvoxreg1,:),1),2));
            mim2(oroi,uroi,:) = squeeze(mean(mean(mim_v(1:nvoxreg1,nvoxreg1+1:end,:),1),2));
            mim2(uroi,oroi,:) = mim2(oroi,uroi,:);        
            mim2(oroi,oroi,:) = squeeze(mean(mean(mim_v(1:nvoxreg1,1:nvoxreg1,:),1),2));
            
        elseif ~strcmp(mode1,'baseline')
            
            %MIC and MIM
            clear a b
            [a , b ] =  fp_mim(Cohroi,npcs([oroi uroi]));
            mic([oroi uroi],[oroi uroi],:) = a;
            mim([oroi uroi],[oroi uroi],:) = b;
        end
    end
end

if strcmp(mode1,'baseline')   
    mic=mic1;
    mim=mim1;
elseif strcmp(mode1,'case2')
    mic = mic2;
    mim=mim2;
elseif strcmp(mode1,'bandc')
    mic.baseline = mic1; 
    mim.baseline = mim1; 
    mic.case2 = mic2; 
    mim.case2 = mim2;    
end

end


