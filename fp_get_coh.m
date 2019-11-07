function coh = fp_get_coh(inname, noEq, vox_ind,flip_id, abs_imag)

load(inname);

if size(
coh(:,:,noEq,:) = [];
match_coh = coh(:,:,vox_ind,:);
flip_coh = match_coh;
flip_coh(:,:,:,4:6) = match_coh(:,:,flip_id,4:6);

%absolute value
if strcmp(abs_imag,'abs')
    r_coh= abs(flip_coh);
elseif strcmp(abs_imag,'imag')
    r_coh = abs(imag(flip_coh));
else
    error('Method unknown!')
end

r_coh(:,1,:,:) = []; %delete inf at freq=1

coh = squeeze(median(r_coh,4));