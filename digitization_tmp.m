nas = [124.4671, 217.2632, 166.6739];
lpa = [ 48.6378, 123.2421, 106.5909];
rpa = [207.9294, 134.2232, 111.0875];

lpa_d =[ 161.0	103.1	125.1];
nas_d = [90.1	172.9	209.5];
rpa_d = [7.0	104.0	119.0];


B = [lpa;nas;rpa]; % points x dimensions
A = [lpa_d;nas_d;rpa_d];

[ret_R, ret_t] = rigid_transform_3D(A, B);

%%

figure; scatter3(A(1,1),A(1,2),A(1,3),'r')
hold on;scatter3(A(2,1),A(2,2),A(2,3),'b')
hold on;scatter3(A(3,1),A(3,2),A(3,3),'g')
hold on; scatter3(B(1,1),B(1,2),B(1,3),'r+')
hold on;scatter3(B(2,1),B(2,2),B(2,3),'b+')
hold on;scatter3(B(3,1),B(3,2),B(3,3),'g+')

hold on; scatter3(B2(1,1),B2(1,2),B2(1,3),'r+')
hold on;scatter3(B2(2,1),B2(2,2),B2(2,3),'b+')
hold on;scatter3(B2(3,1),B2(3,2),B2(3,3),'g+')


%%
u = [Channel.Loc]';
ii = cs_convert(icbm_mri, 'scs', 'mni', u)


%%

A = coord; 

for ii = 1:64
    
    scatter3(A(ii,1),A(ii,2),A(ii,3),'b')
    hold on 
end

A = coord_B; 
for ii = 1:64
    
    scatter3(A(ii,1),A(ii,2),A(ii,3),'r')
    hold on 
end

%%
B2=B_test';
hold on; scatter3(B2(1,1),B2(1,2),B2(1,3),'r+')
hold on;scatter3(B2(2,1),B2(2,2),B2(2,3),'g+')
hold on;scatter3(B2(3,1),B2(3,2),B2(3,3),'b+')

B2=B';
hold on; scatter3(B2(1,1),B2(1,2),B2(1,3),'r.')
hold on;scatter3(B2(2,1),B2(2,2),B2(2,3),'g.')
hold on;scatter3(B2(3,1),B2(3,2),B2(3,3),'b.')

B2 = opoints'; 
hold on; scatter3(B2(1,1),B2(1,2),B2(1,3),'ro')
hold on;scatter3(B2(2,1),B2(2,2),B2(2,3),'go')
hold on;scatter3(B2(3,1),B2(3,2),B2(3,3),'bo')
 
 