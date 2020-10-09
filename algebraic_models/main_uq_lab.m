%% Codes for Sobol's GSA, latent variable approach and Kucherenko approach using UQLab
% simple algebraic models

% see UQLab website: https://www.uqlab.com/

% model 1: Y = X1 + X2 + X2.*X3; flg=1
% model 2: Y = X1 + X2 + X1.*X3; flg=2
% model 3: Y = X1 + X2 + X3 + X4; flg=3
% Xi ~ N(0,1), i = 1:4


clear
close all

addpath(genpath('./'))
% add UQLab path here if necessary

uqlab

flg = 3; % select the model (1,2,3)
rho = 0.7; % define the correlation coefficient between X1 and X4

n_samples = 10000;
param_names = {'X_1','X_2','X_3','X_4','\eta'};
n_params = length(param_names);


%% Parameters distribution

model_fun = @(X) my_reg(X, 0, flg);
model_fun_corr = @(X) my_reg(X, rho, flg); 

ModelOpts.mHandle = model_fun;

ModelOpts.Marginals(1).Name = 'X1';
ModelOpts.Marginals(1).Type = 'Gaussian';
ModelOpts.Marginals(1).Parameters = [0 1];

ModelOpts.Marginals(2).Name = 'X2';
ModelOpts.Marginals(2).Type = 'Gaussian';
ModelOpts.Marginals(2).Parameters = [0 1];

ModelOpts.Marginals(3).Name = 'X3';
ModelOpts.Marginals(3).Type = 'Gaussian';
ModelOpts.Marginals(3).Parameters = [0 1];

ModelOpts.Marginals(4).Name = 'X4';
ModelOpts.Marginals(4).Type = 'Gaussian';
ModelOpts.Marginals(4).Parameters = [0 1];

ModelOpts_nc = ModelOpts;

ModelOpts.Marginals(5).Name = 'eps';
ModelOpts.Marginals(5).Type = 'Gaussian';
ModelOpts.Marginals(5).Parameters = [0 1];


% model dependence for Kucherenko indices
copula = eye(4);
copula(1,4) = rho;
copula(4,1) = rho;

%% Sobol indexes - no correlation

myModel = uq_createModel(ModelOpts_nc);
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

%% Sobol - latent variable

ModelOpts.mHandle = model_fun_corr;
myModel_c = uq_createModel(ModelOpts);
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


%% Kucherenko indices

% setup copula
ModelOpts_nc.Copula.Type = 'Gaussian';
ModelOpts_nc.Copula.Parameters = copula;

myInput = uq_createInput(ModelOpts_nc);
myModel = uq_createModel(ModelOpts_nc);

KucherenkoSensOpts.Type = 'Sensitivity';
KucherenkoSensOpts.Method = 'Kucherenko';
KucherenkoSensOpts.Kucherenko.SampleSize = n_samples;
KucherenkoSensOpts.Model = myModel;

% sensitivity analysis

n_samples_v = [50:50:250, 250:250:5000, 5000:500:10000];
%n_samples_v = 10000;
n_rep = length(n_samples_v);

KucherenkoAnalysis = cell(n_rep,1);
kuch_main = zeros(n_rep,n_params-1);
kuch_total = zeros(n_rep,n_params-1);

for i = 1:n_rep
    KucherenkoSensOpts.Kucherenko.SampleSize = n_samples_v(i);
    KucherenkoAnalysis_i = uq_createAnalysis(KucherenkoSensOpts);
    kuch_main(i,:) = KucherenkoAnalysis_i.Results.FirstOrder';
    kuch_total(i,:) = KucherenkoAnalysis_i.Results.Total';
    KucherenkoAnalysis{i} = KucherenkoAnalysis_i;
end

KuchAnalysis = KucherenkoAnalysis{end};

uq_print(KucherenkoAnalysis{end})
uq_display(KucherenkoAnalysis{end})


%% convergence plots for Kucherenko indices

filename = ['convergence_plot_kucherenko_rho', num2str(rho*10),'_model_',num2str(flg)];

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


%% create output table

filename = ['./results/table_results_rho_', num2str(rho*10),'_model_',num2str(flg),'.xlsx'];

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




