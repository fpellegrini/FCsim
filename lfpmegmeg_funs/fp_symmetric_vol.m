function [sym_pos, noEq] = fp_symmetric_vol(mni_pos)
%kicks out all nodes that do not have an equivalent on the other hemisphere

sym_pos = round(mni_pos);
noEq = []; 

for i = 1:size(sym_pos,1)
    if ~any(sym_pos(:,1)==-sym_pos(i,1) & sym_pos(:,2)==sym_pos(i,2)...
            & sym_pos(:,3)==sym_pos(i,3))
        noEq = [noEq i];        
    end
end
sym_pos(noEq,:)=[];