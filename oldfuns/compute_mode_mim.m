function [mic,mim] = fp_compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA,ihemi)
keyboard
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