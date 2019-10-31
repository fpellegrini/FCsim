function [clu, total] = fp_get_cluster_findclusters(onoff, match_conn,minnbchan)

[clu, total] = findcluster(onoff' ,match_conn, match_conn, minnbchan);
clu = fp_order_clusters(clu',total);