%% simulation of PBPK midazolam model

clear
close all

addpath(genpath('./'))

%% load parameters

% 1: sex     - uniform (threshold 0.5)
% 2: height  - normal
% 3: BMI     - uniform
% 4: MPPGL   - normal
% 5: CYP3A4  - normal
% 6: CYP3A5  - normal
% 7: eta     - normal

param_names = {'sex', 'height', 'BMI', 'MPPGL', 'CYP3A4', 'CYP3A5', 'eta'};

n_samples = 1000;
n_param = length(param_names);
U = rand(n_samples, n_param);

Zd = makedist('Normal', 0, 1);
Z = icdf(Zd, U);

UZ = Z;
UZ(:,1) = U(:,1);
UZ(:,3) = U(:,3);

flg = 1; % 0 not corr | 1 corr
[P_pbpk, P_drug, P_metabolism] = load_physio_param(UZ, flg);

% <1 express CYP3A5 in all the population | <0 do not express CYP3A5 in all the population
cyp3a5_gen = rand(n_samples,1)<1;

%% simulate the model

dose = 5; % [mg]
tspan = [0 24];

n_eq_pbpk = 17;

X0 = zeros(n_eq_pbpk,1);
X0(end) = dose;

t_out_c = cell(n_samples,1);
X_out_c = cell(n_samples,1);
X_out_conc = cell(n_samples,1);
X_liv_c = cell(n_samples,1);

Y = zeros(n_samples, 1);

for i = 1:n_samples
     
    dX_c = @(t,X) midazolam_pbpk(t, X, P_pbpk{i}, P_drug, P_metabolism(i), cyp3a5_gen(i));
    [t_c, X_c] = ode15s(dX_c, tspan, X0);
    
    V_ven = P_pbpk{i}.V_pbpk_system(end);
    
    X_out_c{i} = X_c(:,end);
    X_out_conc{i} = X_c(:,end)/V_ven/P_drug.BP;
    X_liv_c{i} = X_c(:,14);
    t_out_c{i} = t_c;
    
    Y(i) = trapz(t_c, X_c(:,end))./V_ven/P_drug.BP; % plasma
    
end


%% plot results

perc_inf = 2.5;
perc_sup = 97.5;
perc_inf2 = 25;
perc_sup2 = 75;

[ tgrid, Ymin, Ymax, Y_mean, Y_perc_sup, Y_perc_inf, Y_std, Y_median, Y_CV ] = getMinMaxY( t_out_c, X_out_conc, perc_inf, perc_sup );
[ tgrid2, Ymin2, Ymax2, Y_mean2, Y_perc_sup2, Y_perc_inf2, Y_std2, Y_median2, Y_CV2 ] = getMinMaxY( t_out_c, X_out_conc, perc_inf2, perc_sup2 );

t_area = [tgrid; flipud(tgrid)];
t_area2 = [tgrid2; flipud(tgrid2)];

X_area = [Y_perc_sup; flipud(Y_perc_inf)];
X_area2 = [Y_perc_sup2; flipud(Y_perc_inf2)];

face_color = 'c';
face_color2 = 'r';
face_alpha = 0.5;

figure_size = [0 0 0.5 0.8];
format_img = '-depsc'; % -depsc
resolution_img = '-r250';
units_scale = 'normalized';

h = figure();
set(h,'units',units_scale,'outerposition',figure_size)

%subplot(1,2,1)
hold on
fill(t_area, X_area, face_color, 'FaceAlpha', face_alpha, 'LineStyle', 'none')
fill(t_area2, X_area2, face_color2, 'FaceAlpha', face_alpha, 'LineStyle', 'none')
plot(tgrid, Y_median, 'LineWidth', 2)
xlabel('time [h]')
ylabel('conc [mg/L]')
title('midazolam plasma concentration')
ylim([0 1])
xlim([0 24])
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 18)

legend('95% CI','50% CI','median')


printpdf( h, 'plasma_conc', './results', format_img, resolution_img )


%% plot hist & scatt

col0 = [ 0 0.4470 0.7410]; % stdr light blue matlab
col1 = [0 0 1]; % blue
col2 = [1 128/255 0];
col3 = [51/255 51/255 1];
col4 = [153/255 51/255 1];
col5 = [0.85 0.325 0.0980]; %stdr "red" matlab


opts.FigureNumber = 2;
opts.VariableNames = param_names;
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
opts.max_par_fig = 16; % [max_rows, max_columns]
opts.plotonlyscatter = 0;

X = U(:,1:5);

scatter_stats(X,Y,opts)


%% generate boxplot


Ybox = zeros(n_samples, 2);
flg_corr = [0, 1]; % 0 CYP3A5 not expressed | 1 CYP3A5 expressed
cyp3a5_gen_j = 1; % express CYP3A5

for j = 1:2
    
    flg_corr_j = flg_corr(j);
    [P_pbpk, P_drug, P_metabolism] = load_physio_param(UZ, flg_corr_j);
    
    for i = 1:n_samples
        
        dX_c = @(t,X) midazolam_pbpk(t, X, P_pbpk{i}, P_drug, P_metabolism(i), cyp3a5_gen_j);
        [t_c, X_c] = ode15s(dX_c, tspan, X0);
        
        V_ven = P_pbpk{i}.V_pbpk_system(end);
        
        Ybox(i,j) = trapz(t_c, X_c(:,end))./V_ven/P_drug.BP; % plasma
        
    end
    
end


%%% print boxplot

label_v1 = zeros(n_samples,1);
label_v2 = ones(n_samples,1)*0.52;
label_v = categorical([label_v1;label_v2]);
Ybox_v = [Ybox(:,1); Ybox(:,2)];

h1 = figure();
set(h1,'units',units_scale,'outerposition',figure_size)

boxchart(label_v, Ybox_v)
ylabel('AUC [mg/L \cdot h]')
xlabel('\rho_{3A4,3A5}')
title('midazolam plasma AUC')
set(gca, 'fontsize', 18)


printpdf( h1, 'boxplot', './results', format_img, resolution_img )



