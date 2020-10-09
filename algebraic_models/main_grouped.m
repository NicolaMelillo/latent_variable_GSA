%% Script for variance based GSA with groups using the regression models

% model 1: Y = X1 + X2 + X2.*X3; flg=1
% model 2: Y = X1 + X2 + X1.*X3; flg=2
% model 3: Y = X1 + X2 + X3 + X4; flg=3

clear
close all

addpath(genpath('./'))

%% calculate sensitivity indices

flg = 3; % select the model (1,2,3)
rho = 0.9; % define the correlation between X1 and X4

% define the parameters of the analysis
n_param = 3;
n_samples = 10000;
nboot = 1000;
percentiles = [2.5 97.5];
names_param = {'g(X1,X4)', 'X2', 'X3'};

% extract the design matrix for the sobol's GSA
idx_group = 1;
[ Utot ] = getMatrixABCi( n_param, n_samples, idx_group);

% derive the independent parameters
Zd = makedist('Normal', 0, 1);
Z = icdf(Zd, Utot(:,2:end));

% derive the correllated parameters
mu = [0 0];
sigma = [1 rho; rho 1];
Zc = mvnrnd(mu,sigma,n_samples);
Z2 = Zc(Utot(:,1),:);

% define the model inputs and calculate the model outputs
X = [Z2(:,1), Z, Z2(:,2)];
Y = my_reg(X, 0, flg); % here rho=0 because we don't need to use the latent variable

% calculate the variance based indices
input.plot_parameters_names = names_param;
input.nboot = nboot;
input.center_output = 0;
input.n_samples = n_samples;
input.n_param = n_param;
input.percentiles = percentiles;

input.Y = Y;

GSA_reg = SobolSensitivity(input);

%% plot sensitivity indices

input_p.output_SA = GSA_reg;
input_p.param_extracted = names_param;
input_p.ylimit = [-0.1, 1];
input_p.title_bar = 'Y';
input_p.fs = 1;
input_p.print = 0;
input_p.format = '-dpng';
input_p.print_dir = '.images_pub';
input_p.resolution = '-r150';
SobolGraphics( input_p)

%% create output table

filename = ['./results/table_results_rho_', num2str(rho*10),'_model_',num2str(flg),'.xlsx'];

sobol_main = round(GSA_reg{1,1},2);
sobol_main_min = round(GSA_reg{1,3},2);
sobol_main_max = round(GSA_reg{1,4},2);

sobol_tot = round(GSA_reg{2,1},2);
sobol_tot_min = round(GSA_reg{2,3},2);
sobol_tot_max = round(GSA_reg{2,4},2);


T_mean = table(sobol_main, sobol_tot,'RowNames', names_param);
T_CI = table(sobol_main_min, sobol_main_max, sobol_tot_min, sobol_tot_max, 'RowNames', names_param);

writetable(T_mean, filename, 'Sheet', 3, 'WriteRowNames', true)
writetable(T_CI, filename, 'Sheet', 4, 'WriteRowNames', true)

