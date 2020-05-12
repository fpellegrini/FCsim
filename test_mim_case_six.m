
% load('mim_advanced_case_six_results.mat')

for iroi = 2:68
    for jroi = [1:iroi-1]
        clear mc amc imc bmc corr_mc corr_bmc corr_mm corr_bmm
        mc = sum(sum(MIC{iroi,jroi},3),2); 
        [amc, imc] = sort(mc,'descend');
        bmc = zeros(size(mc));
        bmc(imc(1:5))= mc(imc(1:5));
        
        mm = sum(sum(MIM{iroi,jroi},3),2); 
        [amm, imm] = sort(mm,'descend');
        bmm = zeros(size(mm));
        bmm(imm(1:5))= mm(imm(1:5));

        for ii = 1:68
            for jj= 1:68

                corrs_mc(ii,jj) = corr(mc,GT{ii,jj});
                corrs_bmc(ii,jj) = corr(bmc,GT{ii,jj});
                
                corrs_mm(ii,jj) = corr(mm,GT{ii,jj});
                corrs_bmm(ii,jj) = corr(bmm,GT{ii,jj});
            end
        end

%         figure;
%         subplot(2,2,1)
%         hist(corrs_mc(:),100)
%         hold on 
%         plot(mc_gt(iroi,jroi),10,'r+')
%         title(['MC region ' num2str(iroi) ' to region ' num2str(jroi)])
%         subplot(2,2,2)
%         hist(corrs_bmc(:),100)
%         hold on 
%         plot(bmc_gt(iroi,jroi),10,'r+')
%         title(['BMC region ' num2str(iroi) ' to region ' num2str(jroi)])
%         
%         subplot(2,2,3)
%         hist(corrs_mm(:),100)
%         hold on 
%         plot(mm_gt(iroi,jroi),10,'r+')
%         title(['MIM region ' num2str(iroi) ' to region ' num2str(jroi)])
%         subplot(2,2,4)
%         hist(corrs_bmm(:),100)
%         hold on 
%         plot(bmm_gt(iroi,jroi),10,'r+')
%         title(['BMIM region ' num2str(iroi) ' to region ' num2str(jroi)])

        p_mc(iroi,jroi) = sum(corrs_mc(:)>mc_gt(iroi,jroi))/numel(corrs_mc);
        p_bmc(iroi,jroi) = sum(corrs_mc(:)>mc_gt(iroi,jroi))/numel(corrs_mc);
        
        p_mm(iroi,jroi) = sum(corrs_mm(:)>mm_gt(iroi,jroi))/numel(corrs_mm);
        p_bmm(iroi,jroi) = sum(corrs_mm(:)>mm_gt(iroi,jroi))/numel(corrs_mm);
        
    end 
end

outname = './test_mim_results.mat';
save(outname,'p_mc','p_bmc','p_mm','p_bmm','-v7.3')