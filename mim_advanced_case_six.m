
%% sixth case 
clearvars -except mm_gt mc_gt bmm_gt bmc_gt 

for iroi_seed = 1:68
    for iroi_tar = 1:68
        
        fprintf(['seed ' num2str(iroi_seed) ', tar ' num2str(iroi_tar) '\n'])
        tic
        load('./processed_bs/bs_results.mat')
        smooth_cortex = 0.35;
        
        % signal generation
        
        rng(1)
        % number of ROIs in the Desikan-Kiliany Atlas
        nroi = length(cortex.Atlas(3).Scouts);
        
        %set parameters
        n_trials = 200;
        Lepo = 100;
        N = n_trials*Lepo;
        lag = 5;
        fres = 40;
        filtertype= 'd';
        regu=.000001;
        id_trials_1 = 1:n_trials;
        id_trials_2 = 1:n_trials;
        
        %signal generation
        s1 = randn(nroi, N);
        s1(iroi_tar, :) = circshift(s1(iroi_seed, :), lag, 2);
        
        %SNR = 0.5
        s1([iroi_seed iroi_tar],:) = s1([iroi_seed iroi_tar],:)./ norm(s1([iroi_seed iroi_tar],:),'fro');
        s1([(1:iroi_seed-1), (iroi_seed+1):(iroi_tar-1),(iroi_tar+1):end],:) = s1([(1:iroi_seed-1),...
            (iroi_seed+1):(iroi_tar-1),(iroi_tar+1):end],:)./ ...
            norm(s1([(1:iroi_seed-1),(iroi_seed+1):(iroi_tar-1),(iroi_tar+1):end],:),'fro');
        
        %generate ground truth imaginary coherence
        signal_gt = reshape(s1, nroi, Lepo, n_trials);
        CS_gt = fp_tsdata_to_cpsd(signal_gt,fres,'MT',1:nroi, 1:nroi, id_trials_1, id_trials_2);
        CS_gt(:,:,[1 47:end])=[];
        for ifreq = 1: fres
            clear pow
            pow = real(diag(CS_gt(:,:,ifreq)));
            imCoh_gt(:,:,ifreq) = abs(imag(CS_gt(:,:,ifreq)./ sqrt(pow*pow')));
        end
        gt = sum(sum(imCoh_gt,2),3);
        
        % ROI labels
        labels = {cortex.Atlas(3).Scouts.Label};
        %roi inds
        ind_cortex = [];
        sub_ind_cortex = [];
        ind_roi = {};
        sub_ind_roi = {};
        for iROI = 1:nroi
            ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices;
            ind_cortex = cat(1, ind_cortex, ind_roi{iROI});
            [~, ind_roi_cortex{iROI}, ~] = intersect(ind_cortex, ind_roi{iROI}); %index of roi voxels in ind_cortex
            sub_ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices(1);
            sub_ind_cortex = cat(1,sub_ind_cortex, sub_ind_roi{iROI});
            [~,sub_ind_roi_cortex{iROI},~] =  intersect(sub_ind_cortex, sub_ind_roi{iROI});%only one voxel per region
            
        end
        nroi = length(sub_ind_cortex);
        nvox = length(ind_cortex);
        
        L_save = leadfield;
        
        %leadfield for forward model
        L3 = L_save(:, sub_ind_cortex, :);
        for is=1:nroi
            clear L2
            L2 = L3(:,is,:);
            
            %remove radial orientation
            clear u s
            [u, s, v] = svd(squeeze(L2));
            L_forward(:,is,:) = u(:,:)*s(:,1:2);
        end
        
        ni = size(L_forward,3);
        
        p = randn(ni,1);
        p = p/norm(p);
        
        for in = 1:nroi
            L1 = squeeze(L_forward(:,in,:));
            L_mix(:,in) = L1*p;
        end
        
        %project to sensors and add white noise
        for itrial = 1:n_trials
            clear sig whitenoise noise
            sig = L_mix * signal_gt(:,:,itrial);
            sig = sig ./ norm(sig, 'fro');
            
            %add white noise
            whitenoise = randn(size(sig));
            whitenoise = whitenoise ./ norm(whitenoise, 'fro');
            sig = 0.9*sig + 0.1*whitenoise;
            signal_sensor(:,:,itrial) = sig ./ norm(sig, 'fro');
        end
        
        clear L1 L2 L_forward L_mix L3
        
        %% megmeg pipeline start
        %parameters
        
        id_meg_chan = 1:size(signal_sensor,1);
        nmeg = numel(id_meg_chan);
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
        
        mm_gt(iroi_seed,iroi_tar) = corr(mm,gt);
        mc_gt(iroi_seed,iroi_tar) = corr(mc,gt);
        
        bmm_gt(iroi_seed,iroi_tar) = corr(bmm,gt);
        bmc_gt(iroi_seed,iroi_tar) = corr(bmc,gt);
        
        toc
        
    end
end

%%
% 
% imagesc(mic)
% figure
% imagesc(mim)
% figure;
% plot((mc- mean(mc))./std(mc(:)))
% hold on 
% plot((mm - mean(mm))./std(mm(:)))
% legend('mic','mim')
% grid on 

%

% load cm17
% 
% xx = zeros(1,nroi);
% xx([iroi_seed, iroi_tar])=0.2;
% 
% xx1 = gt;
% pos = cortex_highres.Vertices;
% 
% data_in=zeros(1,length(cortex_highres.Curvature));
% allplots_cortex_BS(cortex_highres, data_in, [min(data_in) max(data_in)],...
%     cm17a,'.', smooth_cortex,['ground_thruth_61_45'],  ...
%     {pos(5,:), ...
%     pos(100, :), ...
%     pos(1000, :)});
% clear data_in
% 
% data_in = mm;
% allplots_cortex_BS(cortex, data_in, [min(data_in) max(data_in)],...
%     cm17a,'.', smooth_cortex,['mim_advanced_' num2str(npcs) '_pcs_61_45']);
% 
% data_in = mc;
% allplots_cortex_BS(cortex, data_in, [min(data_in) max(data_in)],...
%     cm17a,'.', smooth_cortex,['mim_advanced_' num2str(npcs) '_pcs_61_45']);
% 
% % close all
