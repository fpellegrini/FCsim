function [TRGC,GC,DIFFGC, to_save] = fp_compute_mode_gc(mode1, D, npcs, dataroi)

ndim = 3; 

%loop over all roi combinations 
for oroi = 1:D.nroi
    for uroi = oroi+1:D.nroi
        
       clear data
       data = cat(2,dataroi{oroi},dataroi{uroi});
       
        if strcmp(mode1,'baseline')||strcmp(mode1,'bandc')
            %preselect "true" voxels (ind sub_ind_roi_region is the
            %index within the roi. Caution with 2 dimensions! 
            
            nvoxreg1 = size(dataroi{oroi},2)/ndim;
            
            
            %indeces of the true voxels 
            clear ouind1 ouind2 baseline inds
            ouind1= cat(2,(D.sub_ind_roi_region{oroi}.*ndim)-2,...
                (D.sub_ind_roi_region{oroi}.*ndim)-1,(D.sub_ind_roi_region{oroi}.*ndim));
            ouind2= cat(2,(D.sub_ind_roi_region{uroi}.*ndim)-2+(nvoxreg1*ndim),...
                (D.sub_ind_roi_region{uroi}.*ndim)-1+(nvoxreg1*ndim),...
                (D.sub_ind_roi_region{uroi}.*ndim)+(nvoxreg1*ndim));
            
            baseline = data(:,[ouind1 ouind2]);    
            iReg = length(D.sub_ind_roi{1});

            %gc across iReg voxels and ndim of one region
            npcs = repmat(iReg*ndim,size(baseline,2)/(ndim*iReg),1);
            inds = fp_npcs2inds(npcs);
            
             %GC and TRGC and DIFFGC
            [trgc, gc, ~, ~] = fp_data2gc(baseline, inds);
            TRGC1(oroi,uroi,:) = trgc(:,1);
            TRGC1(uroi,oroi,:) = trgc(:,2);
            GC1(oroi,uroi,:) = gc(:,1);
            GC1(uroi,oroi,:) = gc(:,2); 
            DIFFGC1(oroi,uroi,:) = TRGC1(oroi,uroi,:)-TRGC1(uroi,oroi,:);

            
        end
        
        if strcmp(mode1,'case2')|| strcmp(mode1,'bandc')
            
            npcs = repmat(ndim,size(data,2)/(ndim),1);
            inds = fp_npcs2inds(npcs);
            [trgc, gc] = fp_data2gc(data, inds);
%             trgc2(oroi,uroi,:) = squeeze(mean(mean(trgc(1:nvoxreg1,nvoxreg1+1:end,:,1),1),2));
%             trgc2(uroi,oroi,:) = squeeze(mean(mean(trgc(1:nvoxreg1,nvoxreg1+1:end,:,1),1),2));
%             
%           
%             %sum up within rois
%             mic2(oroi,uroi,:) = squeeze(mean(mean(mic_v(1:nvoxreg1,nvoxreg1+1:end,:),1),2));
%             mic2(uroi,oroi,:) = mic2(oroi,uroi,:);
%             mic2(oroi,oroi,:) = squeeze(mean(mean(mic_v(1:nvoxreg1,1:nvoxreg1,:),1),2));
%             mim2(oroi,uroi,:) = squeeze(mean(mean(mim_v(1:nvoxreg1,nvoxreg1+1:end,:),1),2));
%             mim2(uroi,oroi,:) = mim2(oroi,uroi,:);        
%             mim2(oroi,oroi,:) = squeeze(mean(mean(mim_v(1:nvoxreg1,1:nvoxreg1,:),1),2));
%             mean_icoh =[];

            
        elseif ~strcmp(mode1,'baseline')
            clear inds trgc gc 
            inds = fp_npcs2inds(npcs);
            
            %GC and TRGC and DIFFGC
            [trgc, gc, ~, ~] = fp_data2gc(data, inds);
            TRGC(oroi,uroi,:) = trgc(:,1);
            TRGC(uroi,oroi,:) = trgc(:,2);
            GC(oroi,uroi,:) = gc(:,1);
            GC(uroi,oroi,:) = gc(:,2); 
            DIFFGC(oroi,uroi,:) = TRGC(oroi,uroi,:)-TRGC(uroi,oroi,:);
        end
    end
end


if strcmp(mode1,'baseline')   
    TRGC=TRGC1;
    GC=GC1;
    DIFFGC = DIFFGC1;
elseif strcmp(mode1,'case2')
    TRGC = TRGC2;
    GC=GC2;
    DIFFGC = DIFFGC2;
elseif strcmp(mode1,'bandc')
    TRGC.baseline = TRGC1;
    clear TRGC1
    GC.baseline = GC1; 
    clear GC1
    DIFFGC.baseline = DIFFGC1; 
    clear DIFFGC1
    TRGC.case2 = TRGC2; 
    clear TRGC2
    GC.case2 = GC2; 
    clear GC2
    DIFFGC.case2 = DIFFGC2; 
    clear DIFFGC2
end
 
to_save.npcs = npcs; 