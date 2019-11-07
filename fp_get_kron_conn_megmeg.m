function kron_conn = fp_get_kron_conn_megmeg(nfreq)

load('roi_conn_new.mat')
roiconn_s = sparse(roi_conn);
freq_conn = fp_get_freq_conn(nfreq);
freq_conn_s = sparse(freq_conn);%nfreq x nfreq
kron_conn = kron(roiconn_s,roiconn_s);
kron_conn = kron(kron_conn,freq_conn_s);
