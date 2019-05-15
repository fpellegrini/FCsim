function freq_conn = fp_get_freq_conn(nfreq)

nfreq = size(coh,2);

freq_conn = zeros(nfreq,nfreq);
for ifreq = 1:nfreq-1
    freq_conn(ifreq,ifreq+1)=1;
    freq_conn(ifreq+1,ifreq)=1;
    freq_conn(ifreq,ifreq) = 1;
end