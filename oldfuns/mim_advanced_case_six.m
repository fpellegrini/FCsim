
%% sixth case 


for iroi_seed = 1:68
    for iroi_tar = 1:68
        clearvars -except mm_gt mc_gt bmm_gt bmc_gt GT MIM MIC iroi_seed iroi_tar
        
        % ROI labels
        [nroi,nvox,ind_cortex,ind_roi_cortex, sub_ind_cortex, roi2vox, leadfield] = fp_get_Desikan;
        
        %signal generation
        fres = 40;
        n_trials = 200;
        [signal_sensor,gt,L_save] = fp_generate_mim_signal(iroi_seed,iroi_tar,fres,n_trials, sub_ind_cortex, leadfield);
        
        
        %% megmeg pipeline start
        %parameters
        id_trials_1 = 1:n_trials;
        id_trials_2 = 1:n_trials;   
        id_meg_chan = 1:size(signal_sensor,1);
        nmeg = numel(id_meg_chan);
        filtertype= 'd';
        regu=.000001;

        CS = fp_tsdata_to_cpsd(signal_sensor,fres,'WELCH',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
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
        
        
        
        %% case 6 specific part 
        
        croi = 1;
        for aroi = 1:nroi
            
            %project to source level
            clear A_ CSv A2
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
        [mic1,mim1]= fp_mim(Cohroi,npcs);
        
        %%
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
        
        GT{iroi_seed,iroi_tar} = gt; 
        MIC{iroi_seed,iroi_tar} = mic1;
        MIM{iroi_seed,iroi_tar} = mim1;
        
        
        toc
        
    end
end

outname = sprintf('./mim_advanced_case_six_results.mat');
save(outname,'mm_gt','mc_gt','bmm_gt','bmc_gt','GT','MIM','MIC','-v7.3')

%%

load cm17
pos = cortex_highres.Vertices;

pos1 = pos(cortex_highres.Atlas(3).Scouts(iroi_seed).Vertices(1,1),:);
pos2 = pos(cortex_highres.Atlas(3).Scouts(iroi_tar).Vertices(1,1),:);

SurfSmoothIterations = ceil(300 * smooth_cortex * length(cortex_highres.Vertices) / 100000);
vc = tess_smooth(cortex_highres.Vertices, 1, SurfSmoothIterations, ...
    tess_vertconn(cortex_highres.Vertices, cortex_highres.Faces), 1);

data_in= gt; %zeros(1,length(cortex_highres.Curvature));
allplots_cortex_BS(cortex_highres1, data_in, [min(data_in) max(data_in)],...
    cm17a,'.', 0,['ground_thruth'],  ...
    {pos1, ...
    pos2});