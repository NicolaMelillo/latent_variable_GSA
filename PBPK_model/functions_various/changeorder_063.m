function [V_changed] = changeorder_063(original_order,changed_order,V_original)
% Change vector order V_original from original_order to changed_order
% original_order and changed_order are structures of length(V_original)x1
% Each cell is a string.
% changed_order must be a permutation of original_order and there must not
% be repeted names.

s1 = size(original_order);

if s1(1)==length(V_original)
    
    lVo = length(V_original);
    
    V_changed = zeros(lVo,1);
    
    for i=1:lVo
        
        posizione = strcmp(original_order{i},changed_order)';
        V_changed(posizione) = V_original(i);
        
    end
    
end
end

