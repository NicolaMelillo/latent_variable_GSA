function [ output, S, T ] = SobolSensitivity( input )
%SOBOLSENSITIVITY compute variance based sensitivity analysis indices


%%% input
% * input: is a structure containing the following fields
%        .Y: output vector, the output is calculated for each row of the
%            matrix obtained from getMatrixABCi.m
%        .n_param: number of parameters
%        .n_samples: number of samples used in getMatrixABCi.m
%        .center_output: option for centering the output during the indices
%                        computation
%        .percentiles: output percentiles from the bootstrap


%%% output
% * output: a cell array in which the first column is the mean, the second the
%           standard deviation, the third is the lower percentile, the fourth is the
%           upper percentile. The first row is the meain effect, the second the total
%           effect.
% * S: n_param X n_boot matrix of the main effects
% * T: n_param X n_boot matrix of the total effects

%% reorganize parameters

n_param = input.n_param;
nboot = input.nboot;
percentiles = input.percentiles;
center_output = input.center_output;
Y = input.Y;
n_samples = input.n_samples;


%% Sensitivity analysis

trig = 0;
if nboot == 0
    trig = 1;
    nboot = 1;
end

S = zeros(n_param,nboot);
T = zeros(n_param,nboot);

for k=1:nboot % bootstrap
    
    fin=0;
    if trig == 0
        boot_idx = randsample(n_samples,n_samples,true);
    else
        boot_idx = 1:n_samples;
    end
    
    ya = Y(fin+1:fin+n_samples,:); fin = fin+n_samples;
    yb = Y(fin+1:fin+n_samples,:); fin = fin+n_samples;
    
    ya = ya(boot_idx);
    yb = yb(boot_idx);
    
    if center_output==1
        ya = ya-mean(ya);
        yb = yb-mean(yb);
    end
    
    yab = [ya; yb];
    f02 = mean(yab)^2;
    Vy = (1/(length(yab))*(sum(yab.^2))-f02);
    
    for i=1:n_param % parameters
        
        yci = Y(fin+1:fin+n_samples); fin = fin+n_samples;
        yci = yci(boot_idx);
        
        % center the output
        if center_output==1
            yci = yci-mean(yci);
        end
        
        % calculate coefficients
        S(i,k) = ( 1/n_samples*(sum(ya.*yci))- 1/n_samples^2*(sum(ya)*sum(yb)) )/Vy; % first effect
        T(i,k) = 1 - (1/n_samples*(sum(yb.*yci))-f02)/Vy; % total effect
    end
    
end

% calculate mean, std and 95/5 percentile of the indexes derived by
% bootstrapping
[ S_mean, S_hi_prctile, S_lo_prctile, S_std  ] = basicstat( S, percentiles(1), percentiles(2) );
[ T_mean, T_hi_prctile, T_lo_prctile, T_std  ] = basicstat( T, percentiles(1), percentiles(2) );

% output formatted as a cell-array
output = {S_mean,S_std,S_lo_prctile,S_hi_prctile;...
    T_mean,T_std,T_lo_prctile,T_hi_prctile};


end



%% References

% (*1) Global Sensitivity Analysis: The Primer
% A. Saltelli, Marco Ratto, Terry Andres, Francesca Campolongo, Jessica Cariboni, Debora Gatelli, Michaela Saisana, Stefano Tarantola
% Wiley 2008
% ISBN: 978-0-470-05997-5
% pag. 164:166
% Check the errata corrige version!





