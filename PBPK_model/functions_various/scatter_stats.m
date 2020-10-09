function [h,n_bins] = scatter_stats(X, Y, opt)

%% Input

% Nicola Melillo, unipv, 10 May 2019

% X: [n_extractions x n_variables] independent variable (e.g. input parameter)
% Y: [n_extractions x 1] dependent variable (e.g. model output)

% opt: structure with fields... All the options should have a default value!
%
%      .VariableNames: cell array containing all variable names, default X1, X2...
%      .Yname: name of the output variable, displaied in the plots
%      .FigureNumber: self explained
%      .n_bins: number of bins used to calculate scatterplot statistics, default sqrt(n_extractions)
%      .percentiles: percentiles to be plotted, in the format [prc_inf1, prc_sup1, prc_inf2, prc_sup2 ...]
%      .facecolor: cell array containing colours of the shaded area corresponding to the percentiles
%      .facealpha: transparency of the area [0-1]
%      .linestyle: stile of the lines corresponding to the percentiles, default 'none'
%      .centraltend: plot central tendency, 1 mean, 2 median, 0 both
%      .color_centraltend: color of the central tendency, cell array {1} mean, {2} median if centraltend==0, otherwise has only size=1.
%      .linewidth: relative to the centraltend
%      .spline_p: p parameter of the cubic splines used to smoooth the percentiles profiles [0,1]
%      .plotscatter: plot the scatterplot
%      .plottails: plot only the extremes pairs of points correspondent to the extreme values of the X
%      .sizepoints: size of the points, default 5
%      .max_par_fig: [max_rows, max_columns] for each plot!
%      .plotonlyscatter: if equal to 1 plot only the scatterplot

n_samples = size(X,1);
n_param = size(X,2);

%% Check options

opt_not_exists = isempty(opt);

% define standard matlab colors...
col0 = [ 0 0.4470 0.7410]; % stdr light blue matlab
col5 = [0.85 0.325 0.0980]; %stdr "red" matlab

% param names
if opt_not_exists || ~isfield(opt, 'VariableNames')  
    param_names = cell(n_param,1);
    for i = 1:n_param
        param_names{i} = ['X',num2str(i)];
    end
else
    param_names = opt.VariableNames;
end

% Y name
if opt_not_exists || ~isfield(opt, 'Yname')  
    Yname = 'Y';
else
    Yname = opt.Yname;
end

% figure number
if opt_not_exists || ~isfield(opt, 'FigureNumber')
    h =  findobj('type','figure');
    n_fig = length(h) + 1;
    clear h;
else
    n_fig = opt.FigureNumber;
end

% number of bins
if opt_not_exists || ~isfield(opt, 'n_bins')
    n_bins = round(sqrt(n_samples));
else
    n_bins = opt.n_bins;
end

% percentiles
if opt_not_exists || ~isfield(opt, 'percentiles')
    percentiles = [2.5 97.5];
    n_percentiles = 1;
else
    percentiles = opt.percentiles;
    n_percentiles = length(percentiles)/2;
end

% facecolor
if opt_not_exists || ~isfield(opt, 'facecolor')
    facecolor = cell(n_percentiles,1);
    for i = 1:n_percentiles
        facecolor{i} = 'r';
    end
else
    if n_percentiles > length(opt.facecolor)
        facecolor = cell(n_percentiles,1);
    for i = 1:n_percentiles
        facecolor{i} = opt.facecolor;
    end
    else
        facecolor = opt.facecolor;
    end
end

% facealpha
if opt_not_exists || ~isfield(opt, 'facealpha')
    facealpha = 0.3;
else
    facealpha = opt.facealpha;
end

% linestyle
if opt_not_exists || ~isfield(opt, 'linestyle')
    linestyle = 'none';
else
    linestyle = opt.linestyle;
end

% centraltend
if opt_not_exists || ~isfield(opt, 'centraltend')
    centraltend = 2; % 1 mean | 2 median | 0 both
else
    centraltend = opt.centraltend;
end

% color_centraltend
if opt_not_exists || ~isfield(opt, 'color_centraltend')
    color_mean = col5;
    color_median = col0;
else
    if centraltend == 0
        color_mean = opt.color_centraltend{1};
        color_median = opt.color_centraltend{2};
    elseif centraltend == 1
        color_mean = opt.color_centraltend{1};
    elseif centraltend == 2
        color_median = opt.color_centraltend{1};
    end
end

% linewidth
if opt_not_exists || ~isfield(opt, 'linewidth')
    linewidth = 2;
else
    linewidth = opt.linewidth;
end

% parameters p and w of the spline...
if opt_not_exists || ~isfield(opt, 'spline_p')
    spline_p = 1;
else
    spline_p = opt.spline_p;
end

% plot scatter
if opt_not_exists || ~isfield(opt, 'plotscatter')
    plotscatter = 0;
else
    plotscatter = opt.plotscatter;
end

% plot tails
if opt_not_exists || ~isfield(opt, 'plottails')
    plottails = 0;
else
    if plotscatter == 1 % no need to plot the tails if the scatter is plotted...
        plottails = 0;
    else
        plottails = opt.plottails;
    end
end

% size point tails
if opt_not_exists || ~isfield(opt, 'sizepoints')
    sizepoints = 5;
else
    sizepoints = opt.sizepoints;
end

% max_subplots
if opt_not_exists || ~isfield(opt, 'max_par_fig')
    max_par_fig = 12;
else
    max_par_fig = opt.max_par_fig;
end

% plot only scatter
if opt_not_exists || ~isfield(opt, 'plotonlyscatter')
    plotonlyscatter = 0;
else
    plotonlyscatter = opt.plotonlyscatter;
    if plotonlyscatter == 1
        plotscatter = 1;
        plottails = 0;
    end
end

% unitscale
if opt_not_exists || ~isfield(opt, 'UnitScale')
    UnitScale = 'centimeters';
else
    UnitScale = opt.UnitScale;
end

% FigureSize
if opt_not_exists || ~isfield(opt, 'FigureSize')
    FigureSize = [0 0 17.4 17.4];
else
    FigureSize = opt.FigureSize;
end

%% Divide all variables into bins...

n_Xbin = floor(n_samples/n_bins);
n_Xleft = n_samples - n_Xbin*n_bins;
size_bins = ones(n_bins,1)*n_Xbin;

for i = 1:n_Xleft % add one data to the first n_Xleft bins...
    size_bins(i) = size_bins(i) + 1;
end

Ybin = cell(n_param,1);
Xbin = cell(n_param,1);

if plottails == 1
    Xij1 = cell(n_param,1);
    Yij1 = cell(n_param,1);
    Xijend = cell(n_param,1);
    Yijend = cell(n_param,1);
end

for i = 1:n_param
    
    Xi = X(:,i);
    [Xis, idx_sort] = sort(Xi);
    Yis = Y(idx_sort);
    
    fin = 0;
    
    Ystats = zeros(n_bins, 2 + length(percentiles));
    Xbin_med = zeros(n_bins,1);
    
    for j = 1:n_bins
        idx_j = fin+1 : fin+size_bins(j);
        Xijs = Xis(idx_j);
        Yijs = Yis(idx_j);
        
        % calculate statistics
        Ystats(j,:) = basicstat_f(Yijs', percentiles);
        Xbin_med(j) = median(Xijs);
        
        fin = fin + size_bins(j);
    end
    
    Xbin{i} = Xbin_med;
    Ybin{i} = Ystats;
    
    % find tails i.e. Xi < median first bin & Xi > median last bin
    % plot all points of the two extreme classes?? look what happens...
    if plottails == 1
        idx_1 = 1 : 1+size_bins(1);
        idx_1 = Xis(idx_1) < median(Xis(idx_1));
        Xij1{i} = Xis(idx_1);
        Yij1{i} = Yis(idx_1);
        
        fin = fin - size_bins(end);
        idx_end = fin : fin+size_bins(end);
        Xis_end = Xis(idx_end);
        Yis_end = Yis(idx_end);
        idx_end = Xis(idx_end) > median(Xis(idx_end));
        Xijend{i} = Xis_end(idx_end);
        Yijend{i} = Yis_end(idx_end);
    end
end

%% plot...

fig_needed = ceil(n_param/(max_par_fig));
n_par_last_fig = mod(n_param, max_par_fig);

n_param_fig = ones(fig_needed,1)*max_par_fig;
if n_par_last_fig ~= 0
    n_param_fig(end) = n_par_last_fig;
end

i2 = 1;
n_fig_k = n_fig-1;

for k = 1:fig_needed
    
    
    [n_col, n_row] = optimizeGrid(n_param_fig(k));
    
    p = spline_p;
    
    h(k) = figure(n_fig_k+1);
    for i = 1:n_param_fig(k)
        Xi = Xbin{i2};
        Yi = Ybin{i2};
        
        pp_mean = csaps(Xi,Yi(:,1),p,Xi);
        pp_median = csaps(Xi,Yi(:,2),p,Xi);
        
        fin = 2;
        pp_inf = cell(n_percentiles,1);
        pp_sup = cell(n_percentiles,1);
        
        subplot(n_row, n_col,i)
        hold on
        
        % plot scatter
        if plotscatter == 1
            scatter(X(:,i2), Y, sizepoints, 'fill');
        end
        
        % plot tails
        if plottails == 1
            scatter(Xij1{i2}, Yij1{i2}, sizepoints, 'fill', 'MarkerFaceColor', col0);
            scatter(Xijend{i2}, Yijend{i2}, sizepoints, 'fill', 'MarkerFacecolor',col0);
        end
        
        if plotonlyscatter == 0 % if I don't want to plot only scatterplot...
            
            % plot percentiles
            for j = 1:n_percentiles
                
                pp_inf{j} = csaps(Xi,Yi(:,fin+1),p,Xi);
                pp_sup{j} = csaps(Xi,Yi(:,fin+2),p,Xi);
                
                Xarea_i = [Xi; flipud(Xi)];
                Yarea_i = [pp_inf{j}; flipud(pp_sup{j})];
                
                fill(Xarea_i, Yarea_i, facecolor{j}, 'FaceAlpha', facealpha, 'LineStyle', linestyle)
                fin = fin+2;
            end
            
            % plot central tendence...
            if centraltend == 2 || centraltend == 0
                plot(Xi, pp_median,'LineWidth',linewidth, 'Color', color_median)
            end
            if centraltend == 1 || centraltend ==0
                plot(Xi, pp_mean,'LineWidth',linewidth, 'Color', color_mean)
            end
            
        end
        
        % labels...
        xlabel(param_names{i2}, 'Interpreter', 'none')
        ylabel(Yname, 'Interpreter', 'none')
        %title_i = [Yname, 32, 'VS', 32, param_names{i2}];
        %title(title_i, 'Interpreter', 'none')
        
        i2 = i2+1;
        
    end
    
    set(h,'units',UnitScale,'outerposition',FigureSize)
    
    n_fig_k = n_fig_k+1;
    
end

end


function [ Sout ] = basicstat_f( S, percentiles)

Smean = mean(S,2);
Smedian = median(S,2);

n_perc = length(percentiles)/2;
Sprc = zeros(1,n_perc*2);
j = 1;
for i = 1:n_perc
    Sprc(j) = prctile(S,percentiles(j),2);
    Sprc(j+1) = prctile(S,percentiles(j+1),2);
    j = j+2;
end

Sout = [Smean, Smedian, Sprc];

end


function [n_cols, n_rows] = optimizeGrid(ls)

% select optimal number of columns & rows for the plots

n_cols = ceil(sqrt(ls));
if (sqrt(ls)-floor(sqrt(ls)))<0.5
    n_rows = floor(sqrt(ls));
else
    n_rows = n_cols;
end

end


