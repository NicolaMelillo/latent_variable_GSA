function [ index ] = findclosest( v_f, v_w )
% Oss: utile per trovare i tempi dei dati (v_f, tipicamente pochi) nelle
%      simulazioni (v_w, tipicamente tante).
% v_f = value to find (data)
% v_w = value where to find (simulation)
% v_f is the vector containg the value to find in v_w
% index is the vector containing the index of the v_f values in v_w, in v_f
% order

l = length(v_f);
index = zeros(l,1);

for i=1:l
   [~,index_i] = min(abs(v_w-v_f(i)));
   index(i) = index_i(1);
end

a = 1;


end

