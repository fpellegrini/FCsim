function fp_megmeg_pipeline_new(patientNumber,DIROUT,filtertype,imethod)

% fp_addpath_sabzi

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

if ~exist('DIROUT','var')
    error('Please indicate where the results should be saved.')
end

DIRLOG = './log/megmeg/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

ndim=2;
nit= 10; %1000
npcs = 3;
fres = 75;

%%
for id = 1:5 %:numel(patientID)
    fprintf('Working on subject %d. \n',id)
    logname = sprintf('%s',patientID{id});
    
    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
        tic
%%        
        %load data
        clear X
        D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
        X = D(:,:,:);
        D_ft = ftraw(D);
        n_trials = length(D_ft.trial);
        
        %channel IDs
        clear id_meg_chan id_lfp_chan
        id_meg_chan = 1:125;
        id_meg_chan(D.badchannels)=[];
        nmeg = numel(id_meg_chan);
        id_lfp_chan = 126:131;
        nlfp = numel(id_lfp_chan);
        
        %scaling
        X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./(10^6);
        
        %construct filters
        
        %load true CS
        clear id_trials_1 id_trials_2 CS
%         id_trials_1 = 1:n_trials;
%         id_trials_2 = 1:n_trials;
%         CS = fp_tsdata_to_cpsd(X,fres,'WELCH',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
    
        load(sprintf('Filter_Patient%s_e.mat',patientID{id}));% 1D-A and true CS
        clear A
        CS = CS(1:(end-nlfp),1:(end-nlfp),:); %throw away lfp channels
        CS(:,:,1)=[];
        
        nfreq = size(CS,3);
        
        %leadfield
        clear L
        load(sprintf('BF_Patient%s.mat',patientID{id}));
        L = fp_get_lf(inverse);
        ns_org = size(L,2);
        
        %get rois
        clear mni_pos label code roi_id u_roi_id csroi
        mni_pos = fp_getMNIpos(patientID{id});
        for ii = 1: ns_org
            [~,~,roi_id(ii)]=fp_get_mni_anatomy_new(mni_pos(ii,:));
        end
        u_roi_id = sort(unique(roi_id));
        nroi = numel(u_roi_id)-1; %because white voxels are not counted 
        
        %construct source filter 
        if strcmp(filtertype,'e')
            A = squeeze(mkfilt_eloreta_v2(L));
            A = permute(A,[1, 3, 2]);
            fqA = ones(1,nfreq);%only one filter for all freqs. 
            nfqA = 1;
            
        elseif strcmp(filtertype,'d')
            A=zeros(nmeg,ndim,ns_org,nfreq);
            
            for ifrq = 1:nfreq
                cCS = CS(:,:,ifrq);
                lambda = mean(diag(real(cCS)))/100;
                
                CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
                
                for is=1:ns_org %iterate across nodes
                    Lloc=squeeze(L(:,is,:));
                    A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
                end
            end
            fqA = 1:nfreq; %This filter is frequency specific. 
            nfqA = nfreq; 
        end
        
        
        clear P
        for aroi = 1:nroi
            
            %project to source level
            clear A_ CSv         
            A_ = A(:, :,roi_id == aroi,:);
            nsroi = size(A_,3);
            A_ = reshape(A_, [nmeg, ndim*nsroi, nfqA]);
            
            for ifq = 1:nfreq
                CSv(ifq,:,:) = A_(:,:,fqA(ifq))' * CS(:,:,ifq) * A_(:,:,fqA(ifq));
            end
            
            %zscoring
            clear ZS CSz
            ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 1))))./diag(squeeze(sum(real(CSv), 1)))));
            for ifreq = 1:nfreq
                CSz(ifreq,:, :) = ZS'*squeeze(CSv(ifreq,:, :))*ZS;
            end
            
            %region pca
            clear CSs v v5
            CSs = squeeze(sum(CSz,1)); %covariance
%             [v, ~, ~] = eig(real(CSs));
%             V{aroi} = v(:,1:npcs); %nregionvoxels*2 x npcs
            
            [V_, D_] = eig(real(CSs));
            [D_, in] = sort(real(diag(D_)), 'descend');
            % variance explained
            vx_ = cumsum(D_)./sum(D_);
            varex{aroi} = vx_;
                      
            
            %concatenate filters
            for ifq = 1:nfqA
                P(:, :, aroi,ifq) = A_(:,:,fqA(ifq)) * ZS * real(V{aroi});
            end
        end
        
        %apply all filters 
        CSroi = [];
        for ifreq = 1:nfreq
            CSroi(:, :, ifreq) = reshape(P, nmeg, [])'*CS(:, :, ifreq)*reshape(P, nmeg, []);
        end
        
        %divide by power to obtain coherence
        clear Cohroi
        for ifreq = 1: nfreq
            clear pow
            pow = real(diag(CSroi(:,:,ifreq)));
            Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
        end
        
        %integrate across npcs
        clear coh
        if strcmp(imethod,'sum')
            
            %sum up coherence across npcs
            ic = 1;
            for iroi = 1:nroi
                jc = 1;
                for jroi = 1:nroi
                    coh(iroi,jroi,:) = squeeze(sum(sum(abs(imag(Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,:))),1),2));
                    jc = jc+npcs;
                end
                ic=ic+npcs;
            end
            
        elseif strcmp(imethod,'mim')
            imethod
        else
            error('Unknown imethod')
        end
           
        
%%      calculate coherence for permutations
       
        for iit = 1:nit
            tic
            fprintf('Working on iteration %d. \n',iit)
            
            clear CS P_shuf
            id_trials_1 = 1:n_trials;
            rng('shuffle')
            id_trials_2 = randperm(n_trials);
            CS = fp_tsdata_to_cpsd(X,fres,'WELCH',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
            
            for aroi = 1: nroi
                
                %project CS to source level
                clear A_ CSv
                A_ = A(:, :,roi_id == aroi,:);
                nsroi = size(A_,3);
                A_ = reshape(A_, [nmeg, ndim*nsroi]);
                
                for ifq = 1:nfreq
                    CSv(ifq,:,:) = A_' * CS(:,:,ifq) * A_;
                end
                
                %zscoring
                clear ZS CSz
                ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 1))))./diag(squeeze(sum(real(CSv), 1)))));
                for ifreq = 1:nfreq
                    CSz(ifreq,:, :) = ZS'*squeeze(CSv(ifreq,:, :))*ZS;
                end
            
                P_shuf(:, :, aroi) = A_*ZS*real(V{aroi});
            end
            
            %apply all filters
            CSroi = [];
            for ifreq = 1:nfreq
                CSroi(:, :, ifreq) = reshape(P_shuf, nmeg, [])'*CS(:, :, ifreq)*reshape(P_shuf, nmeg, []);
            end
            
            %divide by power to obtain coherence
            clear Cohroi
            for ifreq = 1: nfreq
                clear pow
                pow = real(diag(CSroi(:,:,ifreq)));
                Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
            end
            
            %integrate across npcs
            clear coh
            if strcmp(imethod,'sum')
%                 imethod
                
                %sum up coherence across npcs
                ic = 1;
                for iroi = 1:nroi
                    jc = 1;
                    for jroi = 1:nroi
                        coh(iroi,jroi,:) = squeeze(sum(sum(abs(imag(Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,:))),1),2));
                        jc = jc+npcs;
                    end
                    ic=ic+npcs;
                end
                               
            elseif strcmp(imethod,'mim')
                imethod
            else
                error('Unknown imethod')
            end
            
            COH(iit,:,:,:) = coh;
            toc
        end
        
        outname = sprintf('%sroi_coh_sub%s',DIROUT,patientID{id});
        save(outname,'COH','true_coh','-v7.3')
        
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
        toc
    end
end
