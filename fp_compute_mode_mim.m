function [mic,mim,to_save,mean_icoh, mean_acoh] = fp_compute_mode_mim(mode1, D, npcs, V, A2, ZS, CS,fqA,nfqA,ihemi)

ndim = 3; 

%when ihemi == 1: makes sure that rois have the same npcs in both hemispheres
if (strcmp(mode1,'max')||strcmp(mode1,'percent')) && (ihemi==1)    
    npcs = fp_symmetrize_hemispheres(D,npcs,mode1);    
end

%loop over all roi combinations 
for oroi = 1:D.nroi
    for uroi = oroi+1:D.nroi
        
        clear P
        [nmeg, dummy, ~] = size(A2{oroi});
        nvoxreg1 = dummy/ndim;
        nfreq = size(CS,3); 
     
        %concatenate filters
        if strcmp(mode1,'bandc')|| strcmp(mode1,'baseline')||strcmp(mode1,'case2')
            P = cat(2,A2{oroi},A2{uroi});
            P_save{oroi} = A2{oroi}; 
        else
            croi = 1;
            for jroi = [oroi uroi]
                for ifq = 1:nfqA
                    P(:, croi:croi+npcs(jroi)-1,ifq) = A2{jroi}(:,:,fqA(ifq)) * ZS{jroi} * real(V{jroi}(:, 1:npcs(jroi)));
                end
                croi = croi +npcs(jroi);
            end
            P_save{oroi} = P(:,1:npcs(oroi),:);
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
            ouind1= cat(2,(D.sub_ind_roi_region{oroi}.*ndim)-2,...
                (D.sub_ind_roi_region{oroi}.*ndim)-1,(D.sub_ind_roi_region{oroi}.*ndim));
            ouind2= cat(2,(D.sub_ind_roi_region{uroi}.*ndim)-2+(nvoxreg1*ndim),...
                (D.sub_ind_roi_region{uroi}.*ndim)-1+(nvoxreg1*ndim),...
                (D.sub_ind_roi_region{uroi}.*ndim)+(nvoxreg1*ndim));
            
            baseline = Cohroi([ouind1 ouind2],[ouind1 ouind2],:);    
            iReg = length(D.sub_ind_roi{1});

            %mim across iReg voxels and ndim of one region
            npcs = repmat(iReg*ndim,size(baseline,1)/(ndim*iReg),1);
            clear a b
            [a, b]= fp_mim(baseline,npcs);
            mic1([oroi uroi],[oroi uroi],:) = a;
            mim1([oroi uroi],[oroi uroi],:) = b;
            mean_icoh = [];

            
        end
        
        if strcmp(mode1,'case2')|| strcmp(mode1,'bandc')
            
            %MIC and MIM
            fprintf('case2 mim rois %d and %d. \n',oroi,uroi)
            tic
            npcs = repmat(ndim,size(Cohroi,1)/(ndim),1);
            [mic_v,mim_v] =  fp_mim(Cohroi,npcs);
            toc
            
            %sum up mim/ mic within rois
            mic2(oroi,uroi,:) = squeeze(mean(mean(mic_v(1:nvoxreg1,nvoxreg1+1:end,:),1),2));
            mic2(uroi,oroi,:) = mic2(oroi,uroi,:);
            mic2(oroi,oroi,:) = squeeze(mean(mean(mic_v(1:nvoxreg1,1:nvoxreg1,:),1),2));
            mim2(oroi,uroi,:) = squeeze(mean(mean(mim_v(1:nvoxreg1,nvoxreg1+1:end,:),1),2));
            mim2(uroi,oroi,:) = mim2(oroi,uroi,:);        
            mim2(oroi,oroi,:) = squeeze(mean(mean(mim_v(1:nvoxreg1,1:nvoxreg1,:),1),2));
            mean_icoh =[];

            
        elseif ~strcmp(mode1,'baseline')
            
            %MIC and MIM
            clear a b
            [a , b ] =  fp_mim(Cohroi,npcs([oroi uroi]));
            mic([oroi uroi],[oroi uroi],:) = a;
            mim([oroi uroi],[oroi uroi],:) = b;
            
            mean_icoh(oroi,uroi,:) = mean(mean(abs(imag(Cohroi(1:npcs(oroi),npcs(oroi)+1:end,:))),1),2); 
            mean_icoh(uroi,oroi,:) = mean(mean(abs(imag(Cohroi(npcs(oroi)+1:end,1: npcs(oroi)   ,:))),1),2);
            mean_icoh(oroi,oroi,:) = mean(mean(abs(imag(Cohroi(1:npcs(oroi),1:npcs(oroi)    ,:))),1),2);
            
            mean_acoh(oroi,uroi,:) = mean(mean(abs(Cohroi(1:npcs(oroi),npcs(oroi)+1:end,:)),1),2); 
            mean_acoh(uroi,oroi,:) = mean(mean(abs(Cohroi(npcs(oroi)+1:end,1: npcs(oroi)   ,:)),1),2);
            mean_acoh(oroi,oroi,:) = mean(mean(abs(Cohroi(1:npcs(oroi),1:npcs(oroi)    ,:)),1),2);
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
    clear mic1
    mim.baseline = mim1; 
    clear mim1
    mic.case2 = mic2; 
    clear mic2
    mim.case2 = mim2; 
    clear mim2
end

to_save.P = P_save; 
% to_save.Cohroi = Cohroi_save;
% to_save.CSroi = CSroi_save; 
to_save.npcs = npcs; 