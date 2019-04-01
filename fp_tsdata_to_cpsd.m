function S = fp_tsdata_to_cpsd(X,wholeCS,fres,id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials, window,noverlap)

% Estimate cross-power spectral density from time series data.
% If wholeCS is true, the whole cross spectrum will be calculated. If not, only
% the cross spectrum between meg channels and the 6 lfp channels will be
% calculated. 
% Shuffle id_lfp_trials for surrogate images. 

MEG_CHAN = 1:125;
LFP_CHAN = 126:131;
[~,n_times,n_trials] = size(X);
nchan = numel(id_meg_chan) + numel(id_lfp_chan);
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

if wholeCS==1
    S = 0;

    for r = 1:numel(id_meg_trials)
        x_meg = X(:,MEG_CHAN,id_meg_trials(r)); %running trial for meg 
        x_lfp = X(:,LFP_CHAN,id_lfp_trials(r)); %running trial for lfp (potentially permuted)
        x_current = cat(2,x_meg,x_lfp);
        XX = fp_cpsd_welch(x_current,[id_meg_chan id_lfp_chan],[id_meg_chan id_lfp_chan],fres+1,window,noverlap);
        S = S + XX; 
    end
    S = S/n_trials; % average across trials; formerly: S = pi*S/N; 

    % now fill other half of cpsd matrix with complex conjugates
    for i = 1:nchan
        for j = i+1:nchan
            S(j,i,:) = conj(S(i,j,:));
        end
    end
    
else
    S = 0;
    for r = 1:n_trials
        x_meg = X(:,MEG_CHAN,id_meg_trials(r)); %running trial for meg 
        x_lfp = X(:,LFP_CHAN,id_lfp_trials(r)); %running trial for lfp (potentially permuted)
        x_current = cat(2,x_meg,x_lfp);
        XX = fp_cpsd_welch(x_current,id_meg_chan,id_lfp_chan,fres+1,window,noverlap);
        S = S + XX; 

    end
    S = S/n_trials; % average across trials; formerly: S = pi*S/N;
   

end
    

  
    


