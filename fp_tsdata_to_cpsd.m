function S = fp_tsdata_to_cpsd(X,S,fres,id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials, window,noverlap)

% Estimate cross-power spectral density from time series data.
% If S is empty, the whole cross spectrum will be calculated. If not, only
% the cross spectrum between meg channels and the 6 lfp channels will be
% calculated and replaced in S. 
% Shuffle id_lfp_trials for surrogate images. 


[~,m,N] = size(X); %n channels, m time points, N trials 
n = numel(id_meg_chan) + numel(id_lfp_chan);
X = demean(X);
X = permute(X,[2 1 3]); % transpose row, col (works for single-trial data too)

nfft = 2*fres;

if ~exist('window','var')
    window = min(m,nfft); % according to Chronux ... by default Matlab 'pwelch' splits data into 8 overlapping segments
end
assert(window <= m,'window cannot be longer than data');

if ~exist('overlap','var')
    noverlap = round(window/2);
end
assert(noverlap < window,'overlap must be shorter than window');

if isempty(S)
    S = 0;

    for r = 1:numel(id_meg_trials)
        x_meg = X(:,id_meg_chan,id_meg_trials(r)); %running trial for meg 
        x_lfp = X(:,id_lfp_chan,id_lfp_trials(r)); %running trial for lfp (potentially permuted)
        x_current = cat(2,x_meg,x_lfp);
        XX = fp_cpsd_welch(x_current,1:size(x_current,2),1:size(x_current,2),fres+1,window,noverlap);
        S = S + XX; 
    end
    S = S/N; % average across trials; formerly: S = pi*S/N; 

    % now fill other half of cpsd matrix with complex conjugates
    for i = 1:n
        for j = i+1:n
            S(j,i,:) = conj(S(i,j,:));
        end
    end
    
else
    temp_S = 0;
    for r = 1:N
        x_meg = X(:,id_meg_chan,id_meg_trials(r)); %running trial for meg 
        x_lfp = X(:,id_lfp_chan,id_lfp_trials(r)); %running trial for lfp (potentially permuted)
        x_current = cat(2,x_meg,x_lfp);
        XX = fp_cpsd_welch(x_current,id_meg_chan,id_lfp_chan,fres+1,window,noverlap);
        temp_S = temp_S + XX; 

    end
    temp_S = temp_S/N; % average across trials; formerly: S = pi*S/N;
    
    S(id_meg_chan,id_lfp_chan,:) = temp_S;
    S(id_lfp_chan,id_meg_chan,:) = permute(temp_S,[2 1 3]);

end
    

  
    


