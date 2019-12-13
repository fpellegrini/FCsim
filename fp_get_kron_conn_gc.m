function kron_conn = fp_get_kron_conn_gc(nfreq)

freq_conn = fp_get_freq_conn(nfreq);
freq_conn_s = sparse(freq_conn);
conn = fp_find_neighbours('04');
match_conn = conn(voxID{1},voxID{1});
conn_s = sparse(match_conn);
kron_conn = kron(conn_s,freq_conn_s);
%this kron_conn matches a onoff tensor of size (nvox,nfreqs) 
