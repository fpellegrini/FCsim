function fp_plot_slices(data, pos, id,isens,idim, crange)
%data should have be of size 1 x nvox 
%pos are the positions in mni format 
%nit is number of iterations (default: 500)
%outname example: 'test.nii'
% keyboard

if sum(data)==0
    error('data is all zero') 
end
cube = nan([91 109 91]); %destination cube

%mask for all data points outside the brain
mask= wjn_read_nii('/Users/franziskapellegrini/Documents/Master/Masterarbeit/MasterThesis/wjn_toolbox/spm12/toolbox/FieldMap/brainmask.nii');
mask = mask.img;
cube(mask==0) = 0;

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

data = data+100;
cube = reshape(data(in), [91, 109, 91]);
cube(mask==0) = 0;

%scale
Vq = cube;

% keyboard
a = crange;

figone(80,60)
for ii = 10:10:80
    subplot(2,4,ii/10)
    imagesc(squeeze(Vq(ii,:,:))')
    caxis(a)
    set(gca,'YDir','Normal','XDir','Normal')
end
suptitle(['subject ' num2str(id)])
% colorbar
outname = sprintf('sub_%d_sagittal_channel%d_dim%d.png',id,isens,idim);
print(outname,'-dpng');
close all

figone(80,60)
for ii = 15:10:85
    subplot(2,4,(ii-5)/10)
    imagesc(squeeze(Vq(:,ii,:))')
    caxis(a)
    set(gca,'YDir','Normal','XDir','Reverse')
end
suptitle(['subject ' num2str(id)])
% colorbar
outname = sprintf('sub_%d_coronal_channel%d_dim%d.png',id,isens,idim);
print(outname,'-dpng');
close all

figone(80,60)
for ii = 15:10:85
    subplot(2,4,(ii-5)/10)
    imagesc(squeeze(Vq(:,:,ii))')
    caxis(a)
    set(gca,'YDir','Normal','XDir','Reverse')
end
suptitle(['subject ' num2str(id)])
% colorbar
outname = sprintf('sub_%d_transverse_channel%d_dim%d.png',id,isens,idim);
print(outname,'-dpng');
close all


    





