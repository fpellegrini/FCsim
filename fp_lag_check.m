function fp_lag_check

iband = [8 12]; % frequency band of interaction in Hz
band_inds = find(frqs >= iband(1) & frqs <= iband(2)); % indices of interacting frequencies

% filters for band and highpass
[bband, aband] = butter(2, iband/fs*2);

s1 = randn(N, 1);
s1 = filtfilt(bband, aband, s1);

for ilag = 1:80

    s2 = circshift(squeeze(s1), ilag);
    x1 = hilbert(s1); 
    x2 = hilbert(s2); 
    
    x3 = abs(imag(x1)-imag(x2));
    x3_(ilag) = mean(x3); 
    
%     plot(x3) 
%     hold on 

end 

plot(1:80,x3_)
xlabel('lag size in samples') 
ylabel('phase difference')
grid on 
title('signals filtered in alpha band (8-12 Hz)')
%%


fs = 100; % sampling rate
fres = fs; % number of frequency bins (= fres + 1)
Nmin = 3; % length of recording in minutes
N = Nmin*60*fs; % total number of samples
Lepo = 2*fres; % epoch length, should be consistent with fres
n_trials = N/Lepo; % number of epochs
frqs = sfreqs(fres, fs); % freqs in Hz
ilag = 50;
s1_save = randn(N, 1);

for ib = 1:length(frqs)-2
    
    iband = [frqs(ib+1) frqs(ib+2)]; % frequency band of interaction in Hz
    band_inds = find(frqs >= iband(1) & frqs <= iband(2)); % indices of interacting frequencies
    
    % filters for band and highpass
    [bband, aband] = butter(2, iband/fs*2);
    
    
    s1 = filtfilt(bband, aband, s1_save);


    s2 = circshift(squeeze(s1), ilag);
    x1 = hilbert(s1); 
    x2 = hilbert(s2); 
    
    x3 = abs(imag(x1)-imag(x2));
    x3_(ib) = mean(x3); 
    
%     plot(x3) 
%     hold on 

end 

plot(frqs(2:end-1),x3_)
xlabel('Freq of signals') 
ylabel('Phase difference')
grid on 
title('signals delayed with 50 samples')