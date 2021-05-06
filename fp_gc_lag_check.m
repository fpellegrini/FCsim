function fp_gc_lag_check 

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
%lags = 10:10:80;
lags = 1:8;
il=1;
for ilag = lags
    
    s2 = circshift(squeeze(s1), ilag); 
    
    im = 1;
    for morder = 10:10:50

        x1 = cat(2,s1,s2);
        clear x
        x = (0.9*x1 + 0.1*noise)'; 
        
        x = reshape(x,2,Lepo,[]);
        conn = data2sctrgcmim(x, fres, morder, 0,0, [], [], {'GC'},0);
        
        gc(il,im,:,1:2) = squeeze(conn.GC);
        gc(il,im,:,3) = squeeze(conn.GC(:,:,1) - conn.GC(:,:,2));
        
        
        im = im+1; 
    end
    
    il = il+1; 
end 
%%

for igc = 1:3
    figure
    figone(30,50)
    for ilag = 1:numel(lags)
        sp(ilag) = subplot(2,numel(lags)/2,ilag);
        
        plot(frqs, squeeze(gc(ilag, :,:,igc))')
        xlabel('Freqs')
        ylabel('GC')
        title(['Lag ' num2str(lags(ilag))])
        
        rectX = iband;
        rectY = ylim([sp(ilag)]);
%         pch = patch(sp(ilag), rectX([1 2 2 1]) , rectY([1 1 2 2]), 'r', ...
%             'EdgeColor', 'none', 'FaceAlpha', 0.1);
        
        legend('morder 10','morder 20','morder 30','morder 40','morder 50','morder 60','morder 70','morder 80','alpha')

    end
end


