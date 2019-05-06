function S = fp_tsdata_to_cpsd_old(X,fres,ind_1, ind_2, id_meg_trials, id_lfp_trials, window,noverlap)

% Estimate cross-power spectral density from time series data between
% channels ind_1 and channels ind_2 (which are lfp channels).
% Shuffle id_lfp_trials for surrogate images. 

[n_chan,n_times,n_trials] = size(X);
MEG_CHAN = 1:125; %belong to id_meg_trials 
REST_CHAN = 126:n_chan; %lfp + rest channels; belong to id_lfp_trials
X = demean(X);
X = permute(X,[2 1 3]); % transpose row, col (works for single-trial data too)

nfft = 2*fres;

if ~exist('window','var')
    window = min(n_times,nfft); % according to Chronux ... by default Matlab 'pwelch' splits data into 8 overlapping segments
end
assert(window <= n_times,'window cannot be longer than data');

if ~exist('overlap','var')
    noverlap = round(window/2);
end
assert(noverlap < window,'overlap must be shorter than window');

S = 0;

for r = 1:numel(id_meg_trials)
    x_meg = X(:,MEG_CHAN,id_meg_trials(r)); %running trial for meg 
    x_lfp = X(:,REST_CHAN,id_lfp_trials(r)); %running trial for lfp (potentially permuted)
    x_current = cat(2,x_meg,x_lfp);
    XX = fp_cpsd_welch(x_current,ind_1,ind_2,fres+1,window,noverlap);
    S = S + XX; 
end
S = S/n_trials; % average across trials; formerly: S = pi*S/N; 
    

    

  
    


