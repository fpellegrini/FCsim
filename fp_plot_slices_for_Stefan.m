function fp_plot_slices_for_Stefan(data, pos, id,outname1, crange,mask)
%data should have be of size 1 x nvox 
%pos are the positions in mni format 
%nit is number of iterations (default: 500)
%outname example: 'test.nii'
% keyboard

if sum(data)==0
    error('data is all zero') 
end

%from mni to world coordinates
[x, y, z, ~] = mni2orFROMxyz(pos(:,1), pos(:,2), pos(:,3),[],'mni');
x = round(x);
y=round(y);
z=round(z);

data(z<1) = [];
x(z<1)=[];
y(z<1)=[];
z(z<1)=[];

if size(pos,1) ==1906
    load('nn_in')
else
    try 
       out1 = sprintf('./nn_in_Sub%d', id);
       load(out1)
    catch
        [xc yc zc] = ndgrid(1:91, 1:109, 1:91);
        XYZc = [reshape(xc, [], 1), reshape(yc, [], 1), reshape(zc, [], 1)];

        for ii = 1:size(XYZc, 1)
           ii
           d = eucl(XYZc(ii, :), [x y z]);
           [mi in(ii)] = min(d); %in = index vector
        end

        out1 = sprintf('./nn_in_Sub%d', id);
        save(out1,'in')
    end
end

%offset data and crange to get an image distinguishable from background 
data = data+100;
a(1) = crange(1)+97;
a(2) = crange(2) +100;

Vq = reshape(data(in), [91, 109, 91]);
Vq(mask==0) = 0;

%% plotting

figure
% figone(80,60)
for ii = 10:10:80
    subplot(2,4,ii/10)
    imagesc(squeeze(Vq(ii,:,:))')
    caxis(a)
    set(gca,'YDir','Normal','XDir','Normal')
end
suptitle(['subject ' num2str(id)])
outname = sprintf('%s_sagittal.png',outname1);
print(outname,'-dpng');

figure
% figone(80,60)
for ii = 15:10:85
    subplot(2,4,(ii-5)/10)
    imagesc(squeeze(Vq(:,ii,:))')
    caxis(a)
    set(gca,'YDir','Normal','XDir','Reverse')
end
suptitle(['subject ' num2str(id)])
outname = sprintf('%s_coronal.png',outname1);
print(outname,'-dpng');

figure
% figone(80,60)
for ii = 15:10:85
    subplot(2,4,(ii-5)/10)
    imagesc(squeeze(Vq(:,:,ii))')
    caxis(a)
    set(gca,'YDir','Normal','XDir','Reverse')
end
suptitle(['subject ' num2str(id)])
outname = sprintf('%s_transverse.png',outname1);
print(outname,'-dpng');


    





