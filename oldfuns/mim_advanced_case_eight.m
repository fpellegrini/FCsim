
%% case eight


for iroi_seed1 = 1:68
    for iroi_tar1 = 1:68
        clearvars -except mm_gt mc_gt bmm_gt bmc_gt GT MIM MIC iroi_seed1 iroi_tar1 IROI_TAR2 IROI_SEED2
        
        iroi_seed2 = randi(68);
        iroi_tar2 = randi(68);
        IROI_SEED2(iroi_seed1,iroi_tar1) = iroi_seed2;
        IROI_TAR2(iroi_seed1,iroi_tar1)=iroi_tar2;
        % ROI labels
        [nroi,nvox,ind_cortex,ind_roi_cortex, sub_ind_cortex, roi2vox, leadfield] = fp_get_Desikan;
        
        %signal generation
        fres = 40;
        n_trials = 200;
        [signal_sensor,gt,L_save] = fp_generate_mim_signal4(iroi_seed1,iroi_tar1,iroi_seed2, iroi_tar2, fres,n_trials, sub_ind_cortex, leadfield);
        
        
        %% megmeg pipeline stared
        %parameters
        id_trials_1 = 1:n_trials;
        id_trials_2 = 1:n_trials;   
        id_meg_chan = 1:size(signal_sensor,1);
        nmeg = numel(id_meg_chan);
        filtertype= 'd';
        regu=.000001;

        CS = fp_tsdata_to_cpsd(signal_sensor,fres,'MT',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
        CS(:,:,[1 47:end])=[];
        nfreq = size(CS,3);
        
        %leadfield backward model
        L3 = L_save(:, ind_cortex, :);
        for is=1:nvox
            clear L2
            L2 = L3(:,is,:);
            
            %remove radial orientation
            clear u s
            [u, s, v] = svd(squeeze(L2));
            L_backward(:,is,:) = u(:,:)*s(:,1:2);
        end
        ni = size(L_backward,3);
        
        %construct source filter
        if strcmp(filtertype,'e')
            A = squeeze(mkfilt_eloreta_v2(L_backward));
            A = permute(A,[1, 3, 2]);
            fqA = ones(1,nfreq);%only one filter for all freqs.
            nfqA = 1;
            
        elseif strcmp(filtertype,'d')
            
            A=zeros(nmeg,ni,nvox,nfreq);
            
            for ifrq = 1:nfreq
                cCS = CS(:,:,ifrq);
                lambda = mean(diag(real(cCS)))/100;
                
                CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
                
                for is=1:nvox %iterate across nodes
                    Lloc=squeeze(L_backward(:,is,:));
                    A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
                end
            end
            fqA = 1:nfreq; %This filter is frequency specific.
            nfqA = nfreq;
            
            
        elseif strcmp(filtertype,'l')
            
            
        end
        
        
        
        %% 
        
        croi = 1;
        for aroi = 1:nroi
            
            %project to source level
            clear A_ CSv
            A_ = A(:, :,ind_roi_cortex{aroi},:);
            nvoxroi = size(A_,3);
            A2 = reshape(A_, [nmeg, ni*nvoxroi, nfqA]);
            
            
            for ifq = 1: nfreq
                CSv(:,:,ifq) = squeeze(A2(:,:,fqA(ifq)))' * CS(:,:,ifq)...
                    * squeeze(A2(:,:,fqA(ifq)));
            end
            
            %zscoring
            clear ZS CSz
            ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
            for ifreq = 1:nfreq
                CSz(ifreq,:, :) = ZS'*squeeze(CSv(:, :,ifreq))*ZS;
            end
            
            clear CSs v v5 in V_ D_
            CSs = squeeze(sum(CSz,1)); %covariance
            [V_, D_] = eig(real(CSs));
            [D_, in] = sort(real(diag(D_)), 'descend');
            % variance explained
            vx_ = cumsum(D_)./sum(D_);
            invx = 1:min(length(vx_), nmeg);
            npcs(aroi) = min(find(vx_>0.9));
            
            V{aroi} = V_(:,in(1:npcs(aroi))); %nregionvoxels*2 x npcs
            
            %     %concatenate filters
            %     for ifq = 1:nfqA
            %         P(:, :, aroi,ifq) = A2(:,:,fqA(ifq)) * ZS * real(V{aroi});
            %     end
            
            for ifq = 1:nfqA
                P(:, croi:croi+npcs(aroi)-1,ifq) = A2(:,:,fqA(ifq)) * ZS * real(V{aroi});
            end
            croi = croi +npcs(aroi);
        end
        
        %%
        %apply all filters
        CSroi = [];
        for ifreq = 1:nfreq
            CSroi(:, :, ifreq) = reshape(P(:,:,fqA(ifreq)), nmeg, [])'*CS(:, :, ifreq)...
                *reshape(P(:,:,fqA(ifreq)), nmeg, []);
        end
        
        %divide by power to obtain coherence
        clear Cohroi
        for ifreq = 1: fres
            clear pow
            pow = real(diag(CSroi(:,:,ifreq)));
            Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
        end
        
        
        %%
        
        clear mim1 mic1
        ic=1;
        for iroi = 1:nroi
            
            jc=1;
            for jroi = 1:nroi
                
                for ifq = 1:nfqA
                    cs_red=[];
                    cs_red{1} = Cohroi(ic:ic+npcs(iroi)-1,ic:ic+npcs(iroi)-1,ifq);
                    cs_red{2} = Cohroi(ic:ic+npcs(iroi)-1,jc:jc+npcs(jroi)-1,ifq);
                    cs_red{3} = Cohroi(jc:jc+npcs(jroi)-1,jc:jc+npcs(jroi)-1,ifq);
                    
                    caainv=inv(real(cs_red{1})+regu*eye(npcs(iroi))*mean(diag(real(cs_red{1}))));
                    cab=imag(cs_red{2});
                    cbbinv=inv(real(cs_red{3})+regu*eye(npcs(jroi))*mean(diag(real(cs_red{3}))));
                    X=cab*cbbinv*cab';
                    % MIM Ewald Eq. 14
                    mim1(iroi,jroi,ifq)=(trace(caainv*X));
                    caainvsqrt=sqrtm(caainv);
                    Y=caainvsqrt*X*caainvsqrt; %Eq. 23
                    [~,s,~]=svd(Y);
                    % MIC
                    mic1(iroi,jroi,ifq)=sqrt(s(1,1));
                end
                
                jc = jc+npcs(jroi);
            end
            
            ic=ic+npcs(iroi);
        end
        
        %
        mic = sum(mic1,3);
        mim = sum(mim1,3);
        
        mc = sum(mic,2);
        mm = sum(mim,2);
        
        [amm, imm] = sort(mm,'descend');
        bmm = zeros(size(mm));
        bmm(imm(1:5))= mm(imm(1:5));
        
        [amc, imc] = sort(mc,'descend');
        bmc = zeros(size(mc));
        bmc(imc(1:5))= mc(imc(1:5));
        
        mm_gt(iroi_seed1,iroi_tar1) = corr(mm,gt);
        mc_gt(iroi_seed1,iroi_tar1) = corr(mc,gt);
        
        bmm_gt(iroi_seed1,iroi_tar1) = corr(bmm,gt);
        bmc_gt(iroi_seed1,iroi_tar1) = corr(bmc,gt);
        
        GT{iroi_seed1,iroi_tar1} = gt; 
        MIC{iroi_seed1,iroi_tar1} = mic1;
        MIM{iroi_seed1,iroi_tar1} = mim1;
        
        
        toc
        
    end
end

outname = sprintf('./mim_advanced_case_seven_results.mat');
save(outname,'mm_gt','mc_gt','bmm_gt','bmc_gt','GT','MIM','MIC','IROI_SEED2','IROI_TAR2','-v7.3')


