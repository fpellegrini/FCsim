function fp_add_em


DIRIN = './mim_sim4/';
if ~exist(DIRIN);mkdir(DIRIN); end

%%
for iname = 1
    
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
    
    if iname ==  2
        iInt = 1;
    elseif iname>2 && iname<6
        iInt = iname;
    else
        switch iname
            case 6
                iReg = 2;
            case 7
                isnr = 0.5;
            case 8
                isnr = 0.9;
            case 9
                iss=0;
            case 10
                iss = 0.25;
            case 11
                iss = 0.75;
            case 12
                iss = 1;
            case 13
                ilag = 1;
            case 14
                ifilt = 'e';
            case 15
                ifilt = 'c';
            case 16 
                ifilt = 'cr';
        end
    end
    
    %%
    for iit= 1:nit
        tic
%         try
            iit
            clearvars -except iInt iReg isnr iss ilag ifilt iit nit iname DIRIN
            
            inname = sprintf('mim_iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
                ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
            
            load([DIRIN inname '.mat'])
            
            for ipip = 1:numel(MIM)
                
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
                [mrr_mic(ipip), pr_mic(ipip),hk_mic(ipip),em3_mic(ipip)] = fp_mrr_hk(MIC_,iroi_seed,iroi_tar,1);
                [~,~,~,em3_mic_s(ipip)] = fp_mrr_hk(MIC_(ind,ind),iroi_seed,iroi_tar,1);
                
                [mrr_mim(ipip), pr_mim(ipip),hk_mim(ipip),em3_mim(ipip)] = fp_mrr_hk(MIM_,iroi_seed,iroi_tar,1);        
                [~,~,~,em3_mim_s(ipip)] = fp_mrr_hk(MIM_(ind,ind),iroi_seed,iroi_tar,1);
                
                if ipip ~= 10
                    [mrr_aCoh(ipip), pr_aCoh(ipip),hk_aCoh(ipip),em3_aCoh(ipip)] = fp_mrr_hk(aCOH_,iroi_seed,iroi_tar,1);                    
                    [~,~,~,em3_aCoh_s(ipip)] = fp_mrr_hk(aCOH_(ind,ind),iroi_seed,iroi_tar,1);
                    
                    [mrr_iCoh(ipip), pr_iCoh(ipip),hk_iCoh(ipip),em3_iCoh(ipip)] = fp_mrr_hk(iCOH_,iroi_seed,iroi_tar,1);
                    [~,~,~,em3_iCoh_s(ipip)] = fp_mrr_hk(iCOH_(ind,ind),iroi_seed,iroi_tar,1);
                    
                    if ipip ~= 11 && ipip ~= 12
                        %absolute value of gc and only triu is considered. Metric neglects
                        %the direction of the interaction
                        [mrr_absgc(ipip), pr_absgc(ipip),hk_absgc(ipip),em3_absgc(ipip)] ...
                            = fp_mrr_hk(abs(DIFFGC_),iroi_seed,iroi_tar,1);
                        [~,~,~,em3_absgc_s(ipip)] ...
                            = fp_mrr_hk(abs(DIFFGC_(ind,ind)),iroi_seed,iroi_tar,1);
                        
                        %only positive part of gc is submitted and the whole matrix is
                        %considered. Metric that is strongly influenced by the direction of
                        %the effect
                        clear pos_diffgc
                        pos_diffgc = DIFFGC_;
                        pos_diffgc(pos_diffgc< 0) = 0;
                        [mrr_posgc(ipip), pr_posgc(ipip),hk_posgc(ipip),em3_posgc(ipip)]...
                            = fp_mrr_hk(pos_diffgc,iroi_seed,iroi_tar,0);
                        [~,~,~,em3_posgc_s(ipip)]= fp_mrr_hk(pos_diffgc(ind,ind),iroi_seed,iroi_tar,0);
                        
                        %wrong directions 
                        clear pos_diffgc_w
                        pos_diffgc_w = -DIFFGC_;
                        pos_diffgc_w(pos_diffgc_w < 0) = 0;
                        [mrr_posgc_w(ipip), pr_posgc_w(ipip),hk_posgc_w(ipip),em3_posgc_w(ipip)]...
                            = fp_mrr_hk(pos_diffgc_w,iroi_seed,iroi_tar,0);
                        [~,~,~,em3_posgc_w_s(ipip)]= fp_mrr_hk(pos_diffgc_w(ind,ind),iroi_seed,iroi_tar,0);
                    end
                    
                end
            end
            
            %save only evaluation parameters
            outname1 = sprintf('%smrr_%s.mat',DIRIN,params.logname);
            save(outname1,...
                'mrr_mic','pr_mic','hk_mic','em3_mic','em3_mic_s',...
                'mrr_mim','pr_mim','hk_mim','em3_mim','em3_mim_s',...
                'mrr_aCoh','pr_aCoh','hk_aCoh','em3_aCoh','em3_aCoh_s',...
                'mrr_iCoh','pr_iCoh','hk_iCoh','em3_iCoh','em3_iCoh_s',...
                'mrr_absgc','pr_absgc','hk_absgc','em3_absgc','em3_absgc_s',...
                'mrr_posgc','pr_posgc','hk_posgc','em3_posgc','em3_posgc_s',...
                'mrr_posgc_w','pr_posgc_w','hk_posgc_w','em3_posgc_w','em3_posgc_w_s',...
                '-v7.3')
%
%             save(outname1,...
%                 'mrr_mic','pr_mic','hk_mic','em3_mic',...
%                 'mrr_mim','pr_mim','hk_mim','em3_mim',...
%                 'mrr_aCoh','pr_aCoh','hk_aCoh','em3_aCoh',...
%                 'mrr_iCoh','pr_iCoh','hk_iCoh','em3_iCoh',...
%                 'mrr_absgc','pr_absgc','hk_absgc','em3_absgc',...
%                 'mrr_posgc','pr_posgc','hk_posgc','em3_posgc',...
%                 'mrr_posgc_w','pr_posgc_w','hk_posgc_w','em3_posgc_w',...
%                 '-v7.3')

%         end
toc
    end
end
