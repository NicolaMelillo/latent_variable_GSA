%% Latent variable approach and Kucherenko method by varying rho with the algebraic models examples

% see UQLab website: https://www.uqlab.com/

addpath(genpath('./'))
% add UQLab path here if necessary
clear
close all

uqlab

%% parameters definition

param_names = {'X_1','X_2','X_3','X_4','\eta'};

n_param = length(param_names);

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

ModelOpts_k = ModelOpts;

ModelOpts.Marginals(5).Name = 'eta';
ModelOpts.Marginals(5).Type = 'Gaussian';
ModelOpts.Marginals(5).Parameters = [0 1];


%% Sobol GSA with varying rho

% model 1: Y = X1 + X2 + X2.*X3; flg==1
% model 2: Y = X1 + X2 + X1.*X3; flg==2
flg = 2;

rho = -0.9:0.1:0.9;
lr = length(rho);

sobol_p_main = zeros(n_param, lr);
sobol_total = zeros(n_param, lr);
kuch_main = zeros(n_param-1, lr);
kuch_total = zeros(n_param-1, lr);

model_fun = @(X) my_reg(X, 0, flg);
ModelOpts_k.mHandle = model_fun;

n_samples = 10000;

for i = 1:lr
    model_fun_corr = @(X) my_reg(X, rho(i), flg); 

    
    % Sobol GSA
    ModelOpts.mHandle = model_fun_corr;
    myModel = uq_createModel(ModelOpts);
    myInput = uq_createInput(ModelOpts);
    SobolSensOpts.Type = 'Sensitivity';
    SobolSensOpts.Method = 'Sobol';
    SobolSensOpts.Sobol.SampleSize = n_samples;
    SobolSensOpts.Sobol.Estimator = 'homma';
    SobolSensOpts.Model = myModel;
    SobolAnalysis = uq_createAnalysis(SobolSensOpts);
    sobol_p_main(:,i) = SobolAnalysis.Results.FirstOrder;
    sobol_total(:,i) = SobolAnalysis.Results.Total;
    clear('myModel','myInput')
    
    % Kucherenko GSA
    copula = eye(4);
    copula(1,4) = rho(i);
    copula(4,1) = rho(i);
    ModelOpts_k.Copula.Type = 'Gaussian';
    ModelOpts_k.Copula.Parameters = copula;
    myModel_k = uq_createModel(ModelOpts_k);
    myInput_k = uq_createInput(ModelOpts_k);
    KucherenkoSensOpts.Type = 'Sensitivity';
    KucherenkoSensOpts.Method = 'Kucherenko';
    KucherenkoSensOpts.Kucherenko.SampleSize = n_samples;
    KucherenkoSensOpts.Model = myModel_k;
    KucherenkoAnalysis = uq_createAnalysis(KucherenkoSensOpts);
    kuch_main(:,i) = KucherenkoAnalysis.Results.FirstOrder;
    kuch_total(:,i) = KucherenkoAnalysis.Results.Total;
    clear('myModel_k','myInput_k', 'copula')
    
end

%save('all_variables_var_rho_2.mat');

%% plot

idx_to_plot = [1 2 3 4 5];
ltp = length(idx_to_plot);

idx_to_plot2 = [1 2 3 4];
ltp2 = length(idx_to_plot2);

colorOrder = get(gca, 'ColorOrder');

figure_size = [0 0 1 0.8];
figure_size2 = [0 0 17.4 17.4];
figure_size3 = [0 0 12.9 12.9];
figure_size4 = [0 0 8.4 10];

format_img = '-depsc';
resolution_img = '-r250';
units_scale = 'normalized';

filename = ['variable_rho_model_',num2str(flg)];

h = figure();
set(h,'units',units_scale,'outerposition',figure_size)

subplot(1,2,1)
hold on
for i = 1:ltp
    plot(rho, sobol_p_main(idx_to_plot(i),:),'--+','LineWidth',1,'Color',colorOrder(i,:))
    plot(rho, sobol_total(idx_to_plot(i),:),':d','LineWidth',1,'Color',colorOrder(i,:))
end
ylim([-.05 .7])
xlabel('\rho')
ylabel('Sensitivity index')
title('latent variable approach')
%legend('\epsilon_1 main','\epsilon_1 total', '\epsilon_4 main','\epsilon_4 total','\eta main','\eta total','Location','EastOutside')
legend('\epsilon_1 main','\epsilon_1 total', 'X_2 main','X_2 total','X_3 main','X_3 total','\epsilon_4 main','\epsilon_4 total','\eta main','\eta total','Location','EastOutside')

set(gca, 'FontSize', 15)

subplot(1,2,2)
hold on
for i = 1:ltp2
    plot(rho, kuch_main(idx_to_plot2(i),:),'--+','LineWidth',1,'Color',colorOrder(i,:))
    plot(rho, kuch_total(idx_to_plot2(i),:),':d','LineWidth',1,'Color',colorOrder(i,:))
end
ylim([-.05 .7])
xlabel('\rho')
ylabel('Sensitivity index')
title('Kucherenko approach')
%legend('X_1 main','X_1 total', 'X_4 main','X_4 total','\eta main','\eta total','Location','EastOutside')
legend('X_1 main','X_1 total', 'X_2 main','X_2 total','X_3 main','X_3 total','X_4 main','X_4 total','\eta main','\eta total','Location','EastOutside')
set(gca, 'FontSize', 15)

printpdf( h, filename, './results', format_img, resolution_img )



