
DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim5/';

filts = {'e','l'};

for ifilt = 1:2
    for iit = 1: 100
        
        load([DIRIN 'mim_iInt2_iReg1_snr06_iss05_lag2_filt' filts{ifilt} '_p_corr_iter' num2str(iit) '.mat'])
        
        %% add cross-interactions 
        
       iroi_seed = [iroi_seed; iroi_seed]; 
       iroi_tar = [iroi_tar; iroi_tar(2); iroi_tar(1)];
        
        %% Evaluate results with the percentile rank
        
        [pr_mic(ipip)] = fp_pr(MIC{3},iroi_seed,iroi_tar,1);
        [pr_mim(ipip)] = fp_pr(MIM{3},iroi_seed,iroi_tar,1);
        
        if ipip ~= 10
            [pr_aCoh(ipip)] = fp_pr(aCOH{3},iroi_seed,iroi_tar,1);
            [pr_iCoh(ipip)] = fp_pr(iCOH{3},iroi_seed,iroi_tar,1);
            
            if ipip ~= 11 && ipip ~= 12
                %absolute value of gc and only triu is considered. Metric neglects
                %the direction of the interaction
                [pr_abstrgc(ipip)] = fp_pr(abs(DIFFGC{3}),iroi_seed,iroi_tar,1);
                [pr_absgc(ipip)] = fp_pr(abs(GC{3}),iroi_seed,iroi_tar,1);
                
                %only positive part of gc is submitted and the whole matrix is
                %considered. Metric that is strongly influenced by the direction of
                %the effect
                clear pos_difftrgc
                pos_difftrgc = DIFFGC{3};
                pos_difftrgc(pos_difftrgc< 0) = 0;
                [pr_postrgc(ipip)]...
                    = fp_pr(pos_difftrgc,iroi_seed,iroi_tar,0);
                
                clear pos_diffgc
                pos_diffgc = GC{3};
                pos_diffgc(pos_diffgc< 0) = 0;
                [pr_posgc(ipip)]...
                    = fp_pr(pos_diffgc,iroi_seed,iroi_tar,0);
                
            end
        end
        
        fprintf('Saving... \n')
        %save all
        outname = sprintf('%smim_%s.mat',DIRIN,params.logname);
        save(outname,'-v7.3')
        
        %save only evaluation parameters
        outname1 = sprintf('%spr_%s.mat',DIRIN,params.logname);
        save(outname1,...
            'pr_mic',...
            'pr_mim',...
            'pr_aCoh',...
            'pr_iCoh',...
            'pr_absgc',...
            'pr_posgc',...
            'pr_abstrgc',...
            'pr_postrgc',...
            '-v7.3')
        
    end
end