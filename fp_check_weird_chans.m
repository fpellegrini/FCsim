function fp_check_weird_chans(inds,D)
%plots the location of indicated channels
%inds=[2 51 54 89 68 70 84 98 102 103 104 117 121];

D_ft = ftraw(D);
chanpos = D_ft.grad.chanpos(:,[2 1 3]);
loc = mk_sensors_plane(chanpos);

data = D(1:125,:,:);
dat = zeros(size(data));
dat(inds,:,:)=1;

close all
showfield_general(dat(:,4,4),loc)



