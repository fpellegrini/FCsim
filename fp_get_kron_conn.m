function kron_conn = fp_get_kron_conn(patientNumber, vox_ind,nfreq)

conn = fp_find_neighbours(patientNumber);
match_conn = conn(vox_ind,vox_ind);

freq_conn = fp_get_freq_conn(nfreq);
conn_s = sparse(match_conn); %nvox x nvox
freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
kron_conn = kron(conn_s,freq_conn_s); % nvox*nfreq = nkron