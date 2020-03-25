

DIROUT = './';

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[pos, ~] = fp_find_commonvox;
load(sprintf('%sDIFFGC',DIROUT));
nit = 1000;

[nsubs, nvox, nsides, nfreqs] = size(DIFFGC); 

a = nan(nsides,nsubs,3,2); %nsides, nsubs, 3 fbands, start+end

a(1,1,1,1) = 2; 
a(1,1,1,2) = 10; 
a(1,1,2,1) = 10; 
a(1,1,2,2) = 25; 

a(1,2,1,1) = 2; 
a(1,2,1,2) = 20; 
a(1,2,2,1) = 20; 
a(1,2,2,2) = 80; 

a(1,3,1,1) = 2; 
a(1,3,1,2) = 20; 
a(1,3,2,1) = 20; 
a(1,3,2,2) = 40; 

a(1,4,1,1) = 2; 
a(1,4,1,2) = 10; 
a(1,4,2,1) = 12; 
a(1,4,2,2) = 25; 

a(1,5,1,1) = 4; 
a(1,5,1,2) = 40; 
a(1,5,2,1) = 40; 
a(1,5,2,2) = 50; 
a(1,5,3,1) = 50; 
a(1,5,3,2) = 90; 

a(1,6,1,1) = 2; 
a(1,6,1,2) = 40; 
a(1,6,2,1) = 40; 
a(1,6,2,2) = 50; 

a(1,7,1,1) = 2; 
a(1,7,1,2) = 15; 
a(1,7,2,1) = 25; 
a(1,7,2,2) = 50; 

a(1,8,1,1) = 2; 
a(1,8,1,2) = 8; 
a(1,8,2,1) = 6; 
a(1,8,2,2) = 25; 
a(1,8,3,1) = 25; 
a(1,8,3,2) = 60; 

a(1,9,1,1) = 2; 
a(1,9,1,2) = 8; 
a(1,9,2,1) = 10; 
a(1,9,2,2) = 25; 

a(1,10,1,1) = 2; 
a(1,10,1,2) = 5; 
a(1,10,2,1) = 6; 
a(1,10,2,2) = 20; 
a(1,10,3,1) = 22; 
a(1,10,3,2) = 50; 

a(1,11,1,1) = 2; 
a(1,11,1,2) = 12; 
a(1,11,2,1) = 13; 
a(1,11,2,2) = 25; 

%side2 

a(2,1,1,1) = 2; 
a(2,1,1,2) = 5; 
a(2,1,2,1) = 6; 
a(2,1,2,2) = 22; 

a(2,2,1,1) = 2; 
a(2,2,1,2) = 10; 
a(2,2,2,1) = 25; 
a(2,2,2,2) = 60; 

a(2,3,1,1) = 2; 
a(2,3,1,2) = 10; 
a(2,3,2,1) = 12; 
a(2,3,2,2) = 40; 
a(2,3,3,1) = 45; 
a(2,3,3,2) = 90; 

a(2,4,1,1) = 2; 
a(2,4,1,2) = 14; 

a(2,5,1,1) = 2; 
a(2,5,1,2) = 20; 

a(2,6,1,1) = 2; 
a(2,6,1,2) = 35; 
a(2,6,2,1) = 45; 
a(2,6,2,2) = 70; 

a(2,7,1,1) = 2; 
a(2,7,1,2) = 13; 
a(2,7,2,1) = 14; 
a(2,7,2,2) = 30; 

a(2,8,1,1) = 2; 
a(2,8,1,2) = 8; 
a(2,8,2,1) = 10; 
a(2,8,2,2) = 24; 

a(2,9,1,1) = 2; 
a(2,9,1,2) = 5; 
a(2,9,2,1) = 6; 
a(2,9,2,2) = 16; 

a(2,10,1,1) = 2; 
a(2,10,1,2) = 16; 
a(2,10,2,1) = 18; 
a(2,10,2,2) =42; 

a(2,11,1,1) = 2; 
a(2,11,1,2) = 24; 

%%

for iside = 1:2
    for isub = 1:nsubs 
        for fband = 1:3
            
            if ~isnan(a(iside,isub,fband,1))
                
                clear cv start floor
                start = floor(a(iside,isub,fband,1)/2);
                to = floor(a(iside,isub,fband,2)/2);
                cv = squeeze(sum(DIFFGC(isub,:,iside,[start to]),4));   

%                 cv = cv + abs(min(cv(:)));
                outname = sprintf('diffgc_%s_%d_fband_%d-%d.nii',patientID{isub},iside,a(iside,isub,fband,1),a(iside,isub,fband,2));
                fp_data2nii(abs(cv),pos,[],outname)
            end
        end
        
    end
    
end


%%

load('DIFFGC_lcmv.mat')
[nsubs, nvox, nsides, nfreqs] = size(DIFFGC); 
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[pos, ~] = fp_find_commonvox;


for iside = 1:2
    figure
    ii=1;
    for isub = 1:nsubs
        clear cf cv
        
        %         cv = squeeze(sum(DIFFGC(isub,:,iside,:),4));
        %         cv = cv + abs(min(cv(:)));
%         outname = sprintf('diffgc_%s_%d.nii',patientID{isub},iside);
%         fp_data2nii(cv,pos,[],outname)
        
        cf = squeeze(sum(DIFFGC(isub,:,iside,:),2));
        subplot(3,4,ii)
        bar(cf)
        xlabel('freqs')
        xticklabels = 0:5:92;
        xticks = linspace(1,length(cf), numel(xticklabels));
        set(gca,'XTick', xticks,'XTickLabel',xticklabels)
        ylim([-50 50])
        grid on
        
        ii = ii +1;
        
    end
    
%     outname1 = sprintf('diffgc_subplot_%d.png',iside);
%     print(outname1,'-dpng');
%     close all
    
end
