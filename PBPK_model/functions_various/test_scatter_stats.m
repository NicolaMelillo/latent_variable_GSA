%%

clear
close all

%%

col0 = [ 0 0.4470 0.7410]; % stdr light blue matlab
col1 = [0 0 1]; % blue
col2 = [1 128/255 0];
col3 = [51/255 51/255 1];
col4 = [153/255 51/255 1];
col5 = [0.85 0.325 0.0980]; %stdr "red" matlab

X = randn(10000,6);

Y = X(:,1).^2 + 0.5*X(:,2).^3 + 2*sqrt(abs(X(:,3)));% + X(:,3).*X(:,4);

opts.FigureNumber = 4;
%opts.VariableNames = {'x1','x2','x3','x4', 'x5'};
opts.percentiles = [2.5 97.5 25 75];
opts.spline_p = 0.9;
%opts.n_bins = 1000;
opts.facealpha = 0.5;
opts.facecolor = {'c','r'};
%opts.linestyle = '-.';
opts.centraltend = 2; % 1 mean | 2 median | 0 both
%opts.color_centraltend = {col5 col2};
opts.plotscatter = 0;
opts.plottails = 1;
opts.sizepoints = 4;
opts.linewidth = 2;
opts.max_par_fig = 4; % [max_rows, max_columns]
opts.plotonlyscatter = 0;

scatter_stats(X,Y,opts)








