DIRIN = './mim_sim5/c_times/';
DIRFIG = './figures/mimsim_ana/mim_sim5/Manuscript/';
for ip = [1]
    
    clear nInteractions nRegionInts SNR noise_mix nlag filtertype hemisym
    [nInteractions,nRegionInts,SNR,noise_mix,nlag,filtertype, dimred] = fp_get_params(ip);
    o = 1;
    for iInt = nInteractions
        for iReg = nRegionInts
            for isnr = SNR
                for iss = noise_mix
                    for ilag = nlag
                        for ifi = 1:numel(filtertype)
                            
                            clear ifilt
                            ifilt = filtertype{ifi};
                            
                            iit=100;
                            %create logfile for parallelization
                            logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_iter%d'...
                                ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,dimred,iit);
                            
                            load([DIRIN logname '.mat'])
                            
                            for ipip = pips
                                T(ip,ipip,o,:) = ta(ipip+3,:);
                            end
                            
                            o=o+1;
                        end
                    end
                end
            end
        end
    end
    
end
%%
b = median(squeeze(T),2); 
c = b([1:8 21 22 10 9]);

    cl = ...
    [repmat([0.8 0.7 0.6],6,1);
    repmat([0.8 0.4 0.5],2,1);
    [0.4 0.6 0.7];
    [0.5 0.7 0.5];
    [0.4 0.5 0.6];
    [0.8 0.8 0.8]];


figure; 
figone(8, 30)
ba = bar(c, 'facecolor', 'flat');
ba.CData = cl;
grid on 
ylabel('computation times [sec]')
xtitles = {'1PC','2PCs', '3PCs','4PCs','5PCs','6PCs','90%',...
                '99%','Mean+FC','Central','FC+mean','TrueVox'};
xticks=1:12; 
set(gca,'xtick',xticks,'xticklabels',xtitles)

outname = [DIRFIG 'c_times.eps'];
print(outname,'-depsc');

