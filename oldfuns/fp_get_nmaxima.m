function signal_max = fp_get_nmaxima(signal,nmax)

[~, ind] = sort(signal(:),'descend');
signal_max = zeros(size(signal));
signal_max(ind(1:nmax))= signal(ind(1:nmax));