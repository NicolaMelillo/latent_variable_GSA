function [ Ytot ] = getMatrixABCi( n_param, n_samples, idx_group )

%%% Versions
% Nicola Melillo 01/2017
%   - fist version
% Nicola Melillo 15/09/2020 
%   - added grouping

%%% inputs
% n_param: total number of factors (single and grouped)
% n_samples: number of samples
% idx_group: vector containing the indices of the grouped variables

%%% output
% Ytot: n_samples X n_param matrix. 
%       n_param-length(idx_group) columns are sampled from uniform distribution
%       length(idx_group) columns are sampled from a vector of indices = 1:n_samples

%%% HP: 
%   - max(n_group)<n_param, length(n_group)<n_param
%   - For all the factors that are not grouped the output is sampled from a uniform distribution
%   - For all the grouped factors the indices from 1 to n_samples are extracted

% From extraction type specification creates nboot A, B and Ci matrices.

total_samples = (n_samples*2 + n_samples*n_param);
Ytot = zeros(total_samples,n_param);
fin = 0;

Ya = rand(n_samples,n_param*2);

if ~isnan(idx_group)
    for i = 1:length(idx_group)
        Ya(:, idx_group(i)) = randsample(n_samples, n_samples);
        Ya(:,idx_group(i) + n_param) = randsample(n_samples, n_samples);
    end
end

A = Ya(:, 1:n_param);
B = Ya(:, n_param+1:n_param*2);
Ytot(fin+1:fin+n_samples, :) = A; fin = fin + n_samples;
Ytot(fin+1:fin+n_samples, :) = B; fin = fin + n_samples;

for i=1:n_param
    Ci = B;
    Ci(:, i) = A(:, i);
    Ytot(fin+1:fin+n_samples, :) = Ci; fin = fin + n_samples;
end



end

