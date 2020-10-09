%% Sensitivity analysis on midazolam rPBPK by using UQLab

% UQLab website: https://www.uqlab.com/

clear
close all

addpath(genpath('./'))
% here add UQLab path if necessary

uqlab

% dose = 5mg
% route: iv
% output: blood AUC

%% model & parameters specification

flg_gen = 0; % 1 CYP3A5 expressed | 0 CYP3A5 not expressed

% define the model
flg = 0; % 0 no corr | 1 corr
model_fun = @(U) PBPK_uqlab(U, flg, flg_gen);
ModelOpts.mHandle = model_fun;
ModelOpts.isVectorized = true;
myModel = uq_createModel(ModelOpts);

% define the model - latent variable
flg_c = 1; % 0 no corr | 1 corr
model_fun_c = @(U) PBPK_uqlab(U, flg_c, flg_gen);
ModelOpts_nc.mHandle = model_fun_c;
ModelOpts_nc.isVectorized = true;
myModel_c = uq_createModel(ModelOpts_nc);

% define the parameters
param_names = {'sex', 'height', 'BMI', 'MPPGL', 'CYP3A4', 'CYP3A5', 'eta'};

corr_l3A4_l3A5 = 0.5228;

n_samples = 10000;
n_param = length(param_names);

ModelOpts.Marginals(1).Name = 'sex';
ModelOpts.Marginals(1).Type = 'Uniform';
ModelOpts.Marginals(1).Parameters = [0 1];

ModelOpts.Marginals(2).Name = 'height';
ModelOpts.Marginals(2).Type = 'Gaussian';
ModelOpts.Marginals(2).Parameters = [0 1];

ModelOpts.Marginals(3).Name = 'BMI';
ModelOpts.Marginals(3).Type = 'Uniform';
ModelOpts.Marginals(3).Parameters = [0 1];

ModelOpts.Marginals(4).Name = 'MPPGL';
ModelOpts.Marginals(4).Type = 'Gaussian';
ModelOpts.Marginals(4).Parameters = [0 1];

ModelOpts.Marginals(5).Name = 'CYP3A4';
ModelOpts.Marginals(5).Type = 'Gaussian';
ModelOpts.Marginals(5).Parameters = [0 1];

ModelOpts.Marginals(6).Name = 'CYP3A5';
ModelOpts.Marginals(6).Type = 'Gaussian';
ModelOpts.Marginals(6).Parameters = [0 1];

ModelOpts_nc = ModelOpts;

ModelOpts.Marginals(7).Name = 'eta';
ModelOpts.Marginals(7).Type = 'Gaussian';
ModelOpts.Marginals(7).Parameters = [0 1];

% define the copula for the Kucherenko analysis
copula = eye(n_param-1);
copula(5,6) = corr_l3A4_l3A5;
copula(6,5) = corr_l3A4_l3A5;

%% Sobol indexes - no correlation

tic
myInput = uq_createInput(ModelOpts_nc);

SobolSensOpts.Type = 'Sensitivity';
SobolSensOpts.Method = 'Sobol';
SobolSensOpts.Sobol.Estimator = 'homma';
SobolSensOpts.Sobol.SampleSize = n_samples;
SobolSensOpts.Model = myModel;
SobolSensOpts.Bootstrap.Replications = 1000;
SobolSensOpts.Bootstrap.Alpha = 0.05;

SobolAnalysis = uq_createAnalysis(SobolSensOpts);

uq_print(SobolAnalysis)
uq_display(SobolAnalysis)

toc

clear 'myInput'
pause(2)

%% Sobol indices - latent variable approach


tic
myInput = uq_createInput(ModelOpts);

SobolSensOpts_c.Type = 'Sensitivity';
SobolSensOpts_c.Method = 'Sobol';
SobolSensOpts_c.Sobol.Estimator = 'homma';
SobolSensOpts_c.Sobol.SampleSize = n_samples;
SobolSensOpts_c.Model = myModel_c;
SobolSensOpts_c.Bootstrap.Replications = 1000;
SobolSensOpts_c.Bootstrap.Alpha = 0.05;

SobolAnalysis_prec = uq_createAnalysis(SobolSensOpts_c);

uq_print(SobolAnalysis_prec)
uq_display(SobolAnalysis_prec)

toc

clear 'myInput'
pause(2)

%% Kucherenko indices

n_samples_v = [50:50:250, 250:250:5000, 5000:500:10000]; 

n_rep = length(n_samples_v);

tic

ModelOpts_nc.Copula.Type = 'Gaussian';
ModelOpts_nc.Copula.Parameters = copula;

myInput = uq_createInput(ModelOpts_nc);

KucherenkoSensOpts.Type = 'Sensitivity';
KucherenkoSensOpts.Method = 'Kucherenko';
KucherenkoSensOpts.Model = myModel;

KucherenkoAnalysis = cell(n_rep,1);
kuch_main = zeros(n_rep,n_param-1);
kuch_total = zeros(n_rep,n_param-1);

for i = 1:n_rep
    KucherenkoSensOpts.Kucherenko.SampleSize = n_samples_v(i);
    KucherenkoAnalysis_i = uq_createAnalysis(KucherenkoSensOpts);
    kuch_main(i,:) = KucherenkoAnalysis_i.Results.FirstOrder';
    kuch_total(i,:) = KucherenkoAnalysis_i.Results.Total';
    KucherenkoAnalysis{i} = KucherenkoAnalysis_i;
end

toc

uq_print(KucherenkoAnalysis{end})
uq_display(KucherenkoAnalysis{end})

%% convergence plots for Kucherenko indices

filename = 'convergence_plot_kucherenko_flg_0';

font_size = 16;
figure_size = [0 0 1 0.8];
format_img = '-dpng'; % -depsc
resolution_img = '-r250';
units_scale = 'normalized';

h = figure();
set(h,'units',units_scale,'outerposition',figure_size)

subplot(1,2,1)
plot(n_samples_v, kuch_main, 'LineWidth', 2)
xlabel('number of sample points')
ylabel('main effect')
ylim([-0.1, 1])
set(gca, 'fontsize', font_size)

subplot(1,2,2)
plot(n_samples_v, kuch_total, 'LineWidth', 2)
xlabel('number of sample points')
ylabel('total effect')
legend(param_names(1:end-1))
ylim([-0.1, 1])
set(gca, 'fontsize', font_size)

printpdf( h, filename, './results', format_img, resolution_img )


%% plot tables

KuchAnalysis = KucherenkoAnalysis{end};

filename = './results/table_results_PBPK_model_flg_1.xlsx';

sobol_main = [round(SobolAnalysis.Results.Bootstrap.FirstOrder.Mean,2);0];
sobol_main_min = [round(SobolAnalysis.Results.Bootstrap.FirstOrder.CI(:,:,1),2);0];
sobol_main_max = [round(SobolAnalysis.Results.Bootstrap.FirstOrder.CI(:,:,2),2);0];

sobol_tot = [round(SobolAnalysis.Results.Bootstrap.Total.Mean,2);0];
sobol_tot_min = [round(SobolAnalysis.Results.Bootstrap.Total.CI(:,:,1),2);0];
sobol_tot_max = [round(SobolAnalysis.Results.Bootstrap.Total.CI(:,:,2),2);0];

prec_main = round(SobolAnalysis_prec.Results.Bootstrap.FirstOrder.Mean,2);
prec_main_min = round(SobolAnalysis_prec.Results.Bootstrap.FirstOrder.CI(:,:,1),2);
prec_main_max = round(SobolAnalysis_prec.Results.Bootstrap.FirstOrder.CI(:,:,2),2);

prec_tot = round(SobolAnalysis_prec.Results.Bootstrap.Total.Mean,2);
prec_tot_min = round(SobolAnalysis_prec.Results.Bootstrap.Total.CI(:,:,1),2);
prec_tot_max = round(SobolAnalysis_prec.Results.Bootstrap.Total.CI(:,:,2),2);

% for Kucherenko indices we have the convergence plots
kuch_main_t = [round(KuchAnalysis.Results.FirstOrder,2);0];
kuch_tot_t = [round(KuchAnalysis.Results.Total,2);0];


T_mean = table(sobol_main, sobol_tot, kuch_main_t, kuch_tot_t, prec_main, prec_tot,'RowNames', param_names);
T_CI = table(sobol_main_min, sobol_main_max, sobol_tot_min, sobol_tot_max, prec_main_min, prec_main_max, prec_tot_min, prec_tot_max,'RowNames', param_names);

writetable(T_mean, filename, 'Sheet', 1, 'WriteRowNames', true)
writetable(T_CI, filename, 'Sheet', 2, 'WriteRowNames', true)

%save('all_var_flg_1.mat')
