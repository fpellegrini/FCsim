function fp_lag_check

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

%%
fs= 100; 
fres = fs; 
frqs = sfreqs(fres, fs); % freqs in Hz

Nmin = 3; % length of recording in minutes
N = Nmin*60*fs;

iband = [8 12]; % frequency band of interaction in Hz
band_inds = find(frqs >= iband(1) & frqs <= iband(2)); % indices of interacting frequencies

% filters for band and highpass
[bband, aband] = butter(2, iband/fs*2);

s1 = randn(N, 1);
s1 = filtfilt(bband, aband, s1);
Lepo = 2*fres;
noise = randn(N,2);

%%
il = 1;
lags=1:8;
for ilag = lags

    s2 = circshift(squeeze(s1), ilag);
    x1 = cat(2,s1,s2);
    clear x
    x = (0.9*x1 + 0.1*noise)';
    
    x = reshape(x,2,Lepo,[]);
    conn = data2sctrgcmim(x, fres, 20, 0,0, [], [], {'CS'},0);

    a(:,1,il) = angle(conn.CS(:,1,2));
    a(:,2,il) = angle(conn.CS(:,2,1));
    
    il = il+1;

end 

%%
figure
for ilag = 1:8
    sp(ilag) = subplot(2,4,ilag);
    plot(frqs,squeeze(a(:,1,ilag)))
    
    rectX = iband;
    rectY = ylim([sp(ilag)]);
    pch = patch(sp(ilag), rectX([1 2 2 1]) , rectY([1 1 2 2]), 'r', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.1);
    
    title(['lag ' num2str(lags(ilag))])
    xlabel('Freqs')
end

%% same but negative 
figure
for ilag = 1:8
    sp(ilag) = subplot(2,4,ilag);
    plot(frqs,squeeze(a(:,2,ilag)))
    
    rectX = iband;
    rectY = ylim([sp(ilag)]);
    pch = patch(sp(ilag), rectX([1 2 2 1]) , rectY([1 1 2 2]), 'r', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.1);
    
    title(['lag ' num2str(lags(ilag))])
    xlabel('Freqs')
end
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