function S = fp_tsdata_to_cpsd(X,fres,method,ind_1, ind_2, id_trials_1, id_trials_2, window,noverlap,nw,ntapers)

% Estimate cross-power spectral density from time series data between
% channels ind_1 and channels ind_2.
% Shuffle id_lfp_trials for surrogate images. 
keyboard
[n_chan,n_times,n_trials] = size(X);
X = demean(X);
X = permute(X,[2 1 3]); % transpose row, col (works for single-trial data too)
nfft = 2*fres;
ind_pow = intersect(ind_1, ind_2);

if ~exist('window','var')
    window = min(n_times,nfft); % according to Chronux ... by default Matlab 'pwelch' splits data into 8 overlapping segments
end
assert(window <= n_times,'window cannot be longer than data');

if ~exist('overlap','var')
    noverlap = round(window/2);
end
assert(noverlap < window,'overlap must be shorter than window');

if nargin < 3 || isempty(method)
   method = 'MT';    % default is multi-taper
end

if ~isequal(ind_1,unique(ind_1)) || ~isequal(ind_2, unique(ind_2))
    error('ind_1 and ind_2 must be unique')
end

if strcmpi(method,'MT')    
      
    nchunks = floor(((n_times-noverlap)/(window-noverlap))); % FFT chunks per channel

    if nargin < 10 || isempty(nw)
        nw = 3;
    end

    if nargin < 11 || isempty(ntapers)
        ntapers = 2*nw -1;
    end

    tapers   = dpss(window,nw,ntapers,'calc'); % Slepian sequences: tapers is a matrix of size window x ntapers
    taparray = tapers(:,:,ones(1,n_chan));
    
    S = 0;
    for r = 1:n_trials % works for single-trial too
        
        %trial shuffling           
        x_original = X(:,:,id_trials_1(r));
        x_perm = X(:,:,id_trials_2(r));
        
        clear s
        s = fp_cpsd_mt(x_original,x_perm,ind_1, ind_2,fres+1,window,noverlap,nchunks,taparray);
        
        %replace power values by same-trial power values to keep them
        %constant 
        if ~isempty(ind_pow)
            clear pow
            pow = pmtm(x_original(:,ind_pow),nw,nfft,2*nfft);
            for ii = 1:numel(ind_pow)
                clear a b 
                a = find(ind_1 == ind_pow(ii));
                b = find(ind_2 == ind_pow(ii));
                s(:,a,b) = pow(:,ii);
            end
        end
        S = S+s;
        
    end
    S = permute(S,[2 3 1])/n_trials;
    
    
elseif strcmpi(method,'WELCH')
    
    S = 0;

    for r = 1:n_trials
        %trial shuffling
        x_original = X(:,:,id_trials_1(r));
        x_perm = X(:,:,id_trials_2(r));

        XX = fp_cpsd_welch(x_original,x_perm,ind_1,ind_2,fres+1,window,noverlap);
        S = S + XX; 
    end
    S = S/n_trials; % average across trials; formerly: S = pi*S/N; 
    
else
    error('unknown method ''%s''',method);
end
    

    





