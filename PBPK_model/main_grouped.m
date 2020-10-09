%% simulation of rPBPK

clear
close all

addpath(genpath('./'))

%% load parameters

n_samples = 10000;

% 1: sex     - uniform (threshold 0.5)
% 2: height  - normal
% 3: BMI     - uniform
% 4: MPPGL   - normal
% 5: CYP3A4  - normal
% 6: CYP3A5  - normal
% 7: eta     - normal

param_names = {'sex', 'height', 'BMI', 'MPPGL', 'g(3A4,3A5)'};
n_param = length(param_names);

% extract samples for Sobol's analysis
idx_group = 5;
[ U ] = getMatrixABCi( n_param, n_samples, idx_group);
n_samples_tot = size(U,1);

% derive correlated distribution
mu = [0 0];
corr_l3A4_l3A5 = 0.5228;
sigma = [1 corr_l3A4_l3A5; corr_l3A4_l3A5 1];
Zc = mvnrnd(mu, sigma, n_samples_tot);
Z_group = Zc(U(:,idx_group),:);

% reconduct to normal distribution for height, MPPGL and eta
Zd = makedist('Normal', 0, 1);
Z = icdf(Zd, U(:,[2 4]));

% construct input vectors
UZ = [U(:,1), Z(:,1), U(:,3), Z(:,2), Z_group(:,1), Z_group(:,2)];

flg = 0; % 0 not corr | 1 corr (latent variable)
[P_pbpk, P_drug, P_metabolism] = load_physio_param(UZ, flg);

cyp3a5_gen = rand(n_samples_tot,1)<1;


%% simulate the model

dose = 5; % [mg]
tspan = [0 24*7];

n_eq_pbpk = 17;

X0 = zeros(n_eq_pbpk,1);
X0(end) = dose;

t_out_c = cell(n_samples_tot,1);
X_out_c = cell(n_samples_tot,1);
X_out_conc = cell(n_samples_tot,1);
X_liv_c = cell(n_samples_tot,1);

Y = zeros(n_samples_tot, 1);

parfor i = 1:n_samples_tot
     
    dX_c = @(t,X) midazolam_pbpk(t, X, P_pbpk{i}, P_drug, P_metabolism(i), cyp3a5_gen(i));
    [t_c, X_c] = ode15s(dX_c, tspan, X0);
    
    V_ven = P_pbpk{i}.V_pbpk_system(end);
    
    X_out_c{i} = X_c(:,end);
    X_out_conc{i} = X_c(:,end)/V_ven/P_drug.BP;
    X_liv_c{i} = X_c(:,14);
    t_out_c{i} = t_c;
    
    Y(i) = trapz(t_c, X_c(:,end))./V_ven/P_drug.BP; % plasma
    
end

%% Compute & plot sensitivity indices

nboot = 1000;
percentiles = [2.5 97.5];

% calculate the variance based indices
input.plot_parameters_names = param_names;
input.nboot = nboot;
input.center_output = 1;
input.n_samples = n_samples;
input.n_param = n_param;
input.percentiles = percentiles;

input.Y = Y;

GSA_pbpk = SobolSensitivity(input);

% plot sensitivity indices
input_p.output_SA = GSA_pbpk;
input_p.param_extracted = param_names;
input_p.ylimit = [-0.1, 1];
input_p.title_bar = 'plasma AUC';
input_p.fs = 1;
input_p.print = 0;
input_p.format = '-dpng';
input_p.print_dir = './images_pub';
input_p.resolution = '-r150';
SobolGraphics( input_p)

%% create output table

filename = './results/table_results_PBPK_model_flg_1.xlsx';

sobol_main = round(GSA_pbpk{1,1},2);
sobol_main_min = round(GSA_pbpk{1,3},2);
sobol_main_max = round(GSA_pbpk{1,4},2);

sobol_tot = round(GSA_pbpk{2,1},2);
sobol_tot_min = round(GSA_pbpk{2,3},2);
sobol_tot_max = round(GSA_pbpk{2,4},2);


T_mean = table(sobol_main, sobol_tot,'RowNames', param_names);
T_CI = table(sobol_main_min, sobol_main_max, sobol_tot_min, sobol_tot_max, 'RowNames', param_names);

writetable(T_mean, filename, 'Sheet', 3, 'WriteRowNames', true)
writetable(T_CI, filename, 'Sheet', 4, 'WriteRowNames', true)


