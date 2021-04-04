
%% second case 
clear all

for iroi_seed = 1:68 
    for iroi_tar = 1:68
        clearvars -except mm_gt mc_gt bmm_gt bmc_gt GT MIM MIC iroi_seed iroi_tar
        tic
        
        % ROI labels
        [nroi,nvox,ind_cortex,ind_roi_cortex, sub_ind_cortex, roi2vox,leadfield] = fp_get_Desikan;
        
        %signal generation
        fres = 40;
        n_trials = 200;
        [signal_sensor,gt,L_save] = fp_generate_mim_signal(iroi_seed,iroi_tar,fres,n_trials, sub_ind_cortex, leadfield);
        gt = gt(roi2vox);

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
        
        %% case two specific part 

        A2 = reshape(A,nmeg,ni*nvox,nfreq);

        for ifq = 1: nfreq
            CSroi(:,:,ifq) = squeeze(A2(:,:,fqA(ifq)))' * CS(:,:,ifq)...
                * squeeze(A2(:,:,fqA(ifq)));
        end
        clear Cohroi

        %divide by power to obtain coherence
        for ifreq = 1: fres
            clear pow
            pow = real(diag(CSroi(:,:,ifreq)));
            Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
        end


        %% MIM only to aggegate 2 dimensions 
        
        chan = ni; 

        for ivox = 1:nvox 

            ic = ((ivox-1)*2)+1:((ivox-1)*2)+2;
            for jvox = 1:nvox
                jc = ((jvox-1)*2)+1:((jvox-1)*2)+2;

                for ifq = 1:nfqA
                    cs_red=[];
                    cs_red{1} = Cohroi(ic,ic,ifq); %Caa
                    cs_red{2} = Cohroi(ic,jc,ifq); %Cab
                    cs_red{3} = Cohroi(jc,jc,ifq); %Cbb

                    caainv=inv(real(cs_red{1})+regu*eye(chan)*mean(diag(real(cs_red{1}))));
                    cab=imag(cs_red{2});
                    cbbinv=inv(real(cs_red{3})+regu*eye(chan)*mean(diag(real(cs_red{3}))));
                    X=cab*cbbinv*cab';
                    % MIM Ewald Eq. 14
                    mim1(ivox,jvox,ifq)=(trace(caainv*X));
                    caainvsqrt=sqrtm(caainv);
                    Y=caainvsqrt*X*caainvsqrt; %Eq. 23
                    [~,s,~]=svd(Y);
                    % MIC
                    mic1(ivox,jvox,ifq)=sqrt(s(1,1));
                end
            end
        end

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


outname = sprintf('./mim_advanced_case_two_results.mat');
save(outname,'mm_gt','mc_gt','bmm_gt','bmc_gt','-v7.3')    
        

%         imagesc(mic)
%         figure
%         imagesc(mim)
%         figure;
%         plot((mc- mean(mc))./std(mc(:)))
%         hold on 
%         plot((mm - mean(mm))./std(mm(:)))
%         legend('mic','mim')
%         grid on 
% 
%         %%
% 
%         a1 = zeros(size(cortex.Vertices,1),1); 
%         a1(ind_cortex) = mc; 
% 
%         load cm17
%         pos = cortex.Vertices;
% 
%         xx = zeros(size(a1));
%         xx([ind_roi{iroi_seed}; ind_roi{iroi_tar}])=0.2;
% 
%         data_in=xx;
%         allplots_cortex_BS(cortex, data_in, [min(data_in) max(data_in)],...
%             cm17a,'.', smooth_cortex,['ground_thruth_' num2str(iroi_seed) '_' num2str(iroi_tar)]);
%         clear data_in
% 
%         a2 = zeros(size(cortex.Vertices,1),1); 
%         a2(ind_cortex) = mm; 
%         data_in = a2;
%         allplots_cortex_BS(cortex, data_in, [min(data_in) max(data_in)],...
%             cm17a,'.', smooth_cortex,['mim_advanced_case_two_' num2str(iroi_seed) '_' num2str(iroi_tar)]);