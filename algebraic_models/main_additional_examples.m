%% Generic script for Sobol's GSA with any user-defined model

clear
close all

addpath './Sobol_sensitivity'

%% calculate sensitivity indices

% define the parameters of the analysis
n_param = 2;
n_samples = 5000;
nboot = 1000;
percentiles = [2.5 97.5];
names_param = {'g(X1,X2)', 'X3'};

% extract the design matrix for the sobol's GSA
idx_group = 1;
[ Utot ] = getMatrixABCi( n_param, n_samples, idx_group);

% derive the independent parameters
Zd = makedist('Normal', 0, 1);
Z = icdf(Zd, Utot(:,end));

% derive the correllated parameters
mu = [1 1];
rho = 0; % change to 0, 0.9, -0.9 to see the difference in the results
sigma = [1 rho; rho 1];
Zc = mvnrnd(mu, sigma, n_samples);
Zg = Zc(Utot(:,1),:);

% compute output (here you can define each output you want)
Y = Zg(:,1) - Zg(:,2) + Z(:,1); % X1-X2+X3

% compute sensitivity indices
input.plot_parameters_names = names_param;
input.nboot = nboot;
input.center_output = 1;
input.correlation = 0;
input.n_samples = n_samples;
input.n_param = n_param;
input.percentiles = percentiles;

input.Y = Y;

GSA_reg = SobolSensitivity(input);

%% plot sensitivity indices

input_p.output_SA = GSA_reg;
input_p.param_extracted = names_param;
input_p.ylimit = [-0.1, 1];
input_p.title_bar = 'plasma dFdC AUC';
input_p.fs = 1;
input_p.print = 0;
input_p.format = '-dpng';
input_p.print_dir = '../../images_pub';
input_p.resolution = '-r150';
SobolGraphics( input_p)

