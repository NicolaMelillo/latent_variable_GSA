function [ S_mean, S_95prctile, S_05prctile, S_std, S_median, S_CV  ] = basicstat( S, prc_inf, prc_sup, row_names )
% matrix where each row is a different variable and each
% column is a different realization of that variable.

%%%
% Version 25 October 2018
% by Nicola Melillo
% BMS lab - University of Pavia
% Copyright 2018
%%%

S_mean = mean(S,2);
S_05prctile = prctile(S,prc_inf,2);
S_95prctile = prctile(S,prc_sup,2);
S_std = std(S,[],2);
S_median = median(S,2);
S_CV = S_std./S_mean;

if nargin>3
    col_names = {'mean','median','std', 'CV' ,'05prctile','95prctile'};
    T = table(S_mean,S_median,S_std,S_CV,S_05prctile,S_95prctile,'RowNames',row_names,'VariableNames',col_names)
end

end

