

%% signal generation
%parameters

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2;

%leadfield
load(sprintf('BF_Patient%s.mat',patientID{id}));
% L1 = inverse.MEG.L;
% ns = numel(L1);
% for is=1:ns
%     L2(:,:,is) = L1{is}./10^-12;
% end
% L = permute(L2,[1 3 2]);
L = fp_get_lf(inverse);
ni = size(L,3);

for iit = 32
    
    clearvars -except d patientID id seed tar iit L ni delay_ mix time_series
    
    if iit == 1 %baseline
        delay = 30; %in samples, is equal to 100 ms
        inode_seed = 2100;
        inode_tar = 1000;
        inodes = [inode_seed inode_tar];
        isens = 120;        
        
        p = [0.5; 0.5];
        p = p/norm(p);
        for in = 1:2
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'e';
        
    elseif iit > 1 && iit < 6 %varying delays
        
        delay = randi(90);
        inode_seed = 2100;
        inode_tar = 1000;
        inodes = [inode_seed inode_tar];
        isens = 120;
        
        p = [0.5; 0.5];
        p = p/norm(p);
        for in = 1:2
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'e';
        
    elseif iit >5 && iit <12 %varying nodes
        delay = 30;
        inode_seed = randi(size(L,2),1);
        inode_tar = randi(size(L,2),1);
        inodes = [inode_seed inode_tar];
        isens = 120;
        
        p = [0.5; 0.5];
        p = p/norm(p);
        for in = 1:2
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'e';
        
    elseif iit >11 && iit <17 %varying time series
        delay = 30;
        inode_seed = 2100;
        inode_tar = 1000;
        inodes = [inode_seed inode_tar];
        isens = randi(125,1);
        
        p = [0.5; 0.5];
        p = p/norm(p);
        for in = 1:2
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'e';
        
    elseif iit >16 && iit < 21 %varying mixing of leadfield 
        delay = 30;
        inode_seed = 2100;
        inode_tar = 1000;
        inodes = [inode_seed inode_tar];
        isens = randi(125,1);
        
        for in = 1:2
            clear p L1
            p = randn(ni, 1);
            p = p / norm(p);
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'e';
        
    elseif iit == 21
        
        delay = 30; %in samples, is equal to 100 ms
        inode_seed = 2100;
        inode_tar = 1000;
        inodes = [inode_seed inode_tar];
        isens = 120;
        
        p = [0.5; 0.5];
        p = p/norm(p);
        for in = 1:2
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'd';
    elseif iit >21 && iit < 31 %varying nodes, but with dics
        
        delay = 30;
        inode_seed = randi(size(L,2),1);
        inode_tar = randi(size(L,2),1);
        inodes = [inode_seed inode_tar];
        isens = 120;
        
        p = [0.5; 0.5];
        p = p/norm(p);
        for in = 1:2
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'd';
    elseif iit ==31
        
        delay = 0;
        inode_seed = 2100;
        inode_tar = 1000;
        inodes = [inode_seed inode_tar];
        isens = 120;
        
        p = [0.5; 0.5];
        p = p/norm(p);
        for in = 1:2
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'd';
        
    else 
        
        delay = 0;
        inode_seed = 2100;
        inode_tar = 1000;
        inodes = [inode_seed inode_tar];
        isens = 120;
        
        p = [0.5; 0.5];
        p = p/norm(p);
        for in = 1:2
            L1 = squeeze(L(:,inodes(in),:));
            L_mix(:,in) = L1*p;
        end
        filtertype= 'e';
        
        
    end

       
    %load MEG time series and CS
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    X = D(:,:,:);
    id_meg_chan = 1:125;
    X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;
    n_trials = size(X,3);
    
    %correlated nodes
    x1 = squeeze(X(isens,:,:)); %random time meg series
    x2 = x1(delay:end,:); %second time series with time delay
    x1 = x1(1:end-delay+1,:);
    xx(1,:,:)= x1;
    xx(2,:,:) = x2;
    
    %uncorrelated nodes
    % for ii =1:n_rand_nodes
    %     ind = randperm(70);
    %     xx_noise(ii,:,:) = squeeze(X(isens,delay:end,ind));
    % end
    
    n_trials = size(x1,2);
    id_meg_trials = 1:n_trials;
   
    
    %project to sensors and add white noise
    for itrial = 1:n_trials
        clear sig sig1 whitenoise noise
        sig1 = L_mix * xx(:,:,itrial);
        sig1 = sig1 ./ norm(sig1, 'fro');
        %     noise = L_noise * xx_noise(:,:,itrial);
        %     noise = noise ./ norm(noise, 'fro');
        %     sig = alph * sig1 + (1-alph)* noise;
        sig=sig1;
        
        %add white noise
        whitenoise = randn(size(sig));
        whitenoise = whitenoise ./ norm(whitenoise, 'fro');
        sig = 0.9*sig + 0.1*whitenoise;
        signal(:,:,itrial) = sig ./ norm(sig, 'fro');
    end
    
    clear x1 x2 xx whitenoise sig X L1 noise
    
    %% megmeg pipeline start
    %parameters
    
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frqs(frqs>90) = [];
    nfreq = numel(frqs);
    
    
    id_meg_chan = 1:size(signal,1);
    nmeg = numel(id_meg_chan);
    
    id_trials_1 = 1:n_trials;
    id_trials_2 = 1:n_trials;
    CS = fp_tsdata_to_cpsd(signal,fres,'WELCH',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
    CS(:,:,[1 47:end])=[];
    nfreq = size(CS,3);
    
    clear L
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    L = fp_get_lf(inverse);
    ndim = size(L,3);
    
    %construct source filter
    if strcmp(filtertype,'e')
        A = squeeze(mkfilt_eloreta_v2(L));
        A = permute(A,[1, 3, 2]);
        fqA = ones(1,nfreq);%only one filter for all freqs.
        nfqA = 1;
        
    elseif strcmp(filtertype,'d')
        ns_org = size(L,2);
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
    elseif strcmp(filtertype,'l')
        
        
    end
    
    A = permute(A,[1 3 2 4]);
    cCS = CS;
    
    for ifq = 1: nfreq
        for idim =1:ndim
            for ifq = 1:nfreq
                CSv(ifq,idim,:) = squeeze(A(:,:,idim,fqA(ifq)))' * cCS(:,:,ifq)...
                    * squeeze(A(:,inode_seed,idim,fqA(ifq)));
            end
            
            %get voxel power
            pv(idim,:,:) = fp_project_power_2D(CS,squeeze(A(:,:,idim,fqA(ifq))))';
            
        end
    end
    
    %     %zscoring
    %     clear ZS CSz
    %     ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 1))))./diag(squeeze(sum(real(CSv), 1)))));
    %     for ifreq = 1:nfreq
    %         CSz(ifreq,:, :) = ZS'*squeeze(CSv(ifreq,:, :))*ZS;
    %     end
    
    coh = CSv;
    for idim = 1:ndim
        for ifreq = 1:nfreq
            coh(ifreq, idim,:, :) = squeeze(CSv(ifreq, idim,:)) ...
                ./ sqrt(squeeze(pv(idim,ifreq,:))*squeeze(pv(idim,ifreq,inode_seed))');
        end
    end
    
    coh1 = squeeze(sum(abs(imag(coh)),2)); % sum across dims
    
    a = squeeze(sum(coh1,1)); % sum across freqs
    
    %%
    
    outname = sprintf('true_target_%d.nii',iit);
    true = zeros(size(a'));
    true(inode_tar)=1;
    fp_data2nii(true,sources.pos,[],outname,id)    
        
    outname = sprintf('seed_%d.nii',iit);
    true = zeros(size(a'));
    true(inode_seed)=1;
    fp_data2nii(true,sources.pos,[],outname,id)
    
    outname = sprintf('map_megmeg_%d.nii',iit);
    fp_data2nii(a',sources.pos,[],outname,id)
    
    %%
    
    inode_rslt = find(a == max(a));
    
    pos = fp_getMNIpos(patientID{id}); 
    tar_pos = pos(inode_tar,:);
    rslt_pos = pos(inode_rslt,:);
    
    d(iit) = eucl(tar_pos, rslt_pos);
    seed(iit) = inode_seed; 
    tar(iit) = inode_tar; 
    delay_(iit) = delay;
    mix(:,iit) = p; 
    time_series(iit) = isens;
    distribution(:,iit) = a; 
    
    
    
end
save('./location_error_mapmegmeg_dics.mat','d','distribution','seed','tar','delay_','mix','time_series','-v7.3')



