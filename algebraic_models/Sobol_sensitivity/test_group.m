%% my code Sobol's GSA

close all
clear all

n_param = 3;
n_samples = 10000;
nboot = 1000;
percentiles = [2.5 97.5];
names_param = {'X1','X2','X3'};

[ Utot ] = getMatrixABCi( n_param, n_samples, []);
Zd = makedist('Normal', 1, 1);
%Z = icdf(Zd, Utot(:,1:2));
Z = icdf(Zd, Utot);

mu = [0 0];
rho = 0;
sigma = [1 rho; rho 1];
Zc = mvnrnd(mu,sigma,n_samples);
%Z2 = Zc(Utot(:,end),:);

%Y = sqrt(4)*Z(:,1) + sqrt(1)*Z(:,2) + sqrt(3)*Z2(:,1) + sqrt(2)*Z2(:,2);

Y = Z(:,1) + Z(:,2).*Z(:,3);

input.plot_parameters_names = names_param;
input.nboot = nboot;
input.center_output = 1;
input.correlation = 0;
input.n_samples = n_samples;
input.n_param = n_param;
input.percentiles = percentiles;

% plasma AUC
idx_out = 1;
input.Y = Y;
input.title_plot = {'plasma dFdC AUC'};
input.plot_sensitivity_indices = 0;
GSA_AUC_dFdC = SobolSensitivity(input);

input_p.output_SA = GSA_AUC_dFdC;
input_p.param_extracted = names_param;
input_p.ylimit = [-Inf, Inf];
input_p.title_bar = 'plasma dFdC AUC';
input_p.fs = 1;
input_p.print = 0;
input_p.format = '-dpng';
input_p.print_dir = '../../images_pub';
input_p.resolution = '-r150';
SobolGraphics( input_p)


