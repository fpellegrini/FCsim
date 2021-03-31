function P = fp_get_pow(X, fres, id_meg_chan, id_lfp_chan, id_meg_trials, id_lfp_trials, window, noverlap)

[~,n_times,n_trials] = size(X);
X = demean(X);
X = permute(X,[2 1 3]); 
nfft = 2*fres;

if ~exist('window','var')
    window = min(n_times,nfft); % according to Chronux ... by default Matlab 'pwelch' splits data into 8 overlapping segments
end
assert(window <= n_times,'window cannot be longer than data');

if ~exist('overlap','var')
    noverlap = round(window/2);
end
assert(noverlap < window,'overlap must be shorter than window');

P=0;
for itrial = 1:n_trials
    x_meg = X(:,id_meg_chan,id_meg_trials(itrial));
    x_lfp = X(:,id_lfp_chan,id_lfp_trials(itrial));
    x_current = cat(2,x_meg,x_lfp);
    [PP, ~] = pwelch(x_current, window, noverlap, nfft);
    P = P+PP; 
    
end
P = P/n_trials;
P=P';