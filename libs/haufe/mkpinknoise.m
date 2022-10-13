function noise=mkpinknoise(n, m, alpha);

% makes m channels of pink noise
% of length n. Each column is one 
% channel

% randn('state',sum(100*clock))
n1=2*ceil((n-1)/2)+1;
scal=sqrt(1./[1:(n1-3)/2].^alpha');
ff=zeros(n1,m);
ff(2:(n1-1)/2,:)=repmat(scal,1,m);
noise=fft(randn(n1,m));
noise=2*real(ifft(noise.*ff));
noise=noise(1:n,:);

return;


