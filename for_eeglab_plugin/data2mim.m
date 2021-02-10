function [mim,mic] = data2mim(signal_sensor,L,ind_roi_cortex)
%Input: 
% signal_sensor in the form of n_sensor x l_epoch x n_trials
% leadfield L in the form of n_sensors x n_voxels x 3 dimensions
% ind_roi_cortex which is a 1xnrois cell array that provides the indices of
% the voxels per region 
%Output:
% MIM in the form of nroi x nroi x nfreq 
% MIC in the form of nroi x nroi x nfreq

[n_sensors, l_epoch, n_trials] = size(signal_sensor);
ndim = size(L,3);
nroi = numel(ind_roi_cortex);
fres = l_epoch/2;
npcs = 5;

%% cross spectrum

fprintf('Calculating cross spectrum... \n')
CS = tsdata_to_cpsd_fast(signal_sensor,fres,'WELCH');
CS_save = CS;
CS(:,:,1)=[];
nfreq = size(CS,3);

%% lcmv

cCS = sum(CS_save,3);
reg = 0.05*trace(cCS)/length(cCS);
Cr = cCS + reg*eye(size(cCS,1));

[~, A] = lcmv(Cr, L, struct('alpha', 0, 'onedim', 0));
A = permute(A,[1, 3, 2]);


%% PCA

for aroi = 1:nroi
    
    %filter at current roi
    clear A_ CSv
    A_ = A(:, :,ind_roi_cortex{aroi},:);
    nvoxroi(aroi) = size(A_,3); %voxels in the current roi
    A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
    
    %project CS to voxel space
    for ifq = 1:nfreq
        CSv(:,:,ifq) =A2{aroi}' * CS(:,:,ifq) * A2{aroi};
    end
    
    %PCA
    clear CSs v v5 in V_ D_
    CSs = squeeze(sum(real(CSv),3)); %covariance
    [V_, D_] = eig(CSs);
    [D_, in] = sort(real(diag(D_)), 'descend');
    V{aroi} = V_(:,in(1:npcs));
    
    %concatenate filters
    P{aroi} = A2{aroi} * real(V{aroi});
    
end


%% calculate MIM/MIC

%loop over all roi combinations
for oroi = 1:nroi
    for uroi = oroi+1:nroi
        
        %apply all filters
        clear Proi
        CSroi = [];
        Proi = cat(2,P{oroi},P{uroi});
        for ifreq = 1:nfreq
            CSroi(:, :, ifreq) = reshape(Proi, n_sensors, [])'*CS(:, :, ifreq)...
                *reshape(Proi, n_sensors, []);
        end
        
        %divide by power to obtain coherence
        clear Cohroi
        for ifreq = 1: nfreq
            clear pow
            pow = real(diag(CSroi(:,:,ifreq)));
            Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
        end
        clear CSroi
        
        %MIC and MIM
        clear a b
        [a , b ] =  roi_mim(Cohroi,[npcs npcs]);
        mic([oroi uroi],[oroi uroi],:) = a;
        mim([oroi uroi],[oroi uroi],:) = b;
        
    end
end
