function fp_add_em_short

% fp_addpath

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim4/';
if ~exist(DIRIN);mkdir(DIRIN); end

%%
for iname = [8]
    
    iname
    
    clearvars -except iname DIRIN
    
    %default paramenters
    nit = 100;
    iInt = 2;
    iReg=1;
    isnr=0.7;
    iss = 0.5;
    ilag=2;
    ifilt='l';
    dimred = 'p';
    
    if iname ==  2
        iInt = 1;
    elseif iname>2 && iname<7
        iInt = iname;
    else
        switch iname
            case 7
                iReg = 2;
            case 8
                isnr = 0.5;
            case 9
                isnr = 0.9;
            case 10
                iss=0;
            case 11
                iss = 0.25;
            case 12
                iss = 0.75;
            case 13
                iss = 1;
            case 14
                ilag = 1;
            case 15
                ifilt = 'e';
            case 16
                ifilt = 'c';
            case 17
                ifilt = 'cr';
            case 18
                ifilt = 'd';
            case 19 
                dimred = 's';
             
        end
    end
    
    %%
    for iit= 1:nit
        tic
        try
            iit
            clearvars -except iInt iReg isnr iss ilag ifilt iit nit iname DIRIN dimred
            
            if iname == 7 || iname == 18 || iname == 19
                inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d_%s'...
                ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit,dimred);
            
            else
                inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
                    ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
            end
            
            load([DIRIN inname '.mat'])
            
            if iname ~= 16
                pips = 1:numel(MIM);
            else 
                pips = [1:3 8];
            end
            
            for ipip = pips
                
                MIC_ = MIC{ipip};
                MIM_ = MIM{ipip};
                aCOH_ = aCOH{ipip};
                iCOH_ = iCOH{ipip};
                if ipip ~= 10 && ipip ~= 11 && ipip ~= 12  && ipip < 21
                    DIFFGC_ = DIFFGC{ipip}; 
                end
                
                ind = randperm(68,68);
                
                % pr is now not area under the precision-recall curve
                % anymore but the new performance metric 
                [pr_mic(ipip)] = fp_mrr_hk_short(MIC_,iroi_seed,iroi_tar,1);
                
                [pr_mim(ipip)] = fp_mrr_hk_short(MIM_,iroi_seed,iroi_tar,1);        
                
                if ipip ~= 10
                    [ pr_aCoh(ipip)] = fp_mrr_hk_short(aCOH_,iroi_seed,iroi_tar,1);                  
                    
                    [pr_iCoh(ipip)] = fp_mrr_hk_short(iCOH_,iroi_seed,iroi_tar,1);
                    
                    if ipip ~= 11 && ipip ~= 12
                        %absolute value of gc and only triu is considered. Metric neglects
                        %the direction of the interaction
                        [ pr_absgc(ipip)] ...
                            = fp_mrr_hk_short(abs(DIFFGC_),iroi_seed,iroi_tar,1);
                        
                        %only positive part of gc is submitted and the whole matrix is
                        %considered. Metric that is strongly influenced by the direction of
                        %the effect
                        clear pos_diffgc
                        pos_diffgc = DIFFGC_;
                        pos_diffgc(pos_diffgc< 0) = 0;
                        [pr_posgc(ipip)]...
                            = fp_mrr_hk_short(pos_diffgc,iroi_seed,iroi_tar,0);
                        
                        %wrong directions 
                        clear pos_diffgc_w
                        pos_diffgc_w = -DIFFGC_;
                        pos_diffgc_w(pos_diffgc_w < 0) = 0;
                        [pr_posgc_w(ipip)]...
                            = fp_mrr_hk_short(pos_diffgc_w,iroi_seed,iroi_tar,0);
                    end
                    
                end
            end
            
            %% asr            
            
%             CS = fp_tsdata_to_cpsd(signal_sensor,fres,'WELCH',1:size(signal_sensor,1),...
%                 1:size(signal_sensor,1), 1:size(signal_sensor,3), 1:size(signal_sensor,3));
% 
%             ASR = log(norm(imag(CS(:)))/norm(real(CS(:))));
            ASR =[];
            
            
            %%
            
            %save only evaluation parameters and ASR 
            outname1 = sprintf('%smrr_%s_pr.mat',DIRIN,params.logname);
%             save(outname1,...
%                 'mrr_mic','pr_mic','hk_mic','em3_mic','em3_mic_s',...
%                 'mrr_mim','pr_mim','hk_mim','em3_mim','em3_mim_s',...
%                 'mrr_aCoh','pr_aCoh','hk_aCoh','em3_aCoh','em3_aCoh_s',...
%                 'mrr_iCoh','pr_iCoh','hk_iCoh','em3_iCoh','em3_iCoh_s',...
%                 'mrr_absgc','pr_absgc','hk_absgc','em3_absgc','em3_absgc_s',...
%                 'mrr_posgc','pr_posgc','hk_posgc','em3_posgc','em3_posgc_s',...
%                 'mrr_posgc_w','pr_posgc_w','hk_posgc_w','em3_posgc_w','em3_posgc_w_s',...
%                 'ASR', '-v7.3')

            if iname == 7 || iname == 18
                mrr_mic = [];
                hk_mic = []; 
                em3_mic=[]; 
                mrr_mim = [];
                hk_mim = []; 
                em3_mim=[]; 
                mrr_iCoh = [];
                hk_iCoh = []; 
                em3_iCoh=[]; 
                mrr_aCoh = [];
                hk_aCoh = []; 
                em3_aCoh=[]; 
                mrr_absgc = [];
                hk_absgc = []; 
                em3_absgc=[];
                mrr_posgc = [];
                hk_posgc = []; 
                em3_posgc=[];
            end

            save(outname1,...
                'mrr_mic','pr_mic','hk_mic','em3_mic',...
                'mrr_mim','pr_mim','hk_mim','em3_mim',...
                'mrr_aCoh','pr_aCoh','hk_aCoh','em3_aCoh',...
                'mrr_iCoh','pr_iCoh','hk_iCoh','em3_iCoh',...
                'mrr_absgc','pr_absgc','hk_absgc','em3_absgc',...
                'mrr_posgc','pr_posgc','hk_posgc','em3_posgc',...
                'ASR','-v7.3')
            toc
    end

    end
end

