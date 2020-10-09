function [parameters_pbpk2, parameters_drug, P_metabolism] = load_physio_param(UZ, flag_corr)

%%%
% Version 15 April 2020
% by Nicola Melillo, UoM
% whole body PBPK model
%%%

n_samples = size(UZ,1);

% volumes [L]
% time [h]
% drug mass [mg]

% U
% 1: sex     - uniform (threshold 0.5)
% 2: height  - normal
% 3: BMI     - uniform
% 4: MPPGL   - normal
% 5: CYP3A4  - normal
% 6: CYP3A5  - normal
% 7: eta     - normal

%% generate a population

% get organs weight fraction and composition and fraction of cardiac output
% directed to the organ
[compNames, fractions, par_cell_FM ] = getFractions();

% generate organs weight
[ SubjParamOutput, pH ] = getPhysiologicalParameters(compNames, fractions, par_cell_FM, UZ(:,1:3));


%% drug related parameters

parameters_drug.pKa = 6.15; % (*9)
parameters_drug.BP = 0.66; % (*10)
parameters_drug.fup = 0.0303; % (*10)
parameters_drug.mw = 325.77; % [g/mol] % (*9)
parameters_drug.logPow = 3.13; % (*9)
parameters_drug.S = 0;
parameters_drug.pH_S = 0;
parameters_drug.PappAB = 0; % [cm/s]
parameters_drug.PappBA = 0;
parameters_drug.type = '?';

parameters_drug.fut = 1/(1+(1-parameters_drug.fup)/parameters_drug.fup*0.5); % keep if not known, otherwhise specify it (Isoherranen)
parameters_drug.logDvow = 1.115*parameters_drug.logPow-1.35; % keep if not known, otherwhise specify it
parameters_drug.Y_hh = 1+10^(parameters_drug.pKa - pH.plasma);
parameters_drug.Pow = 10^parameters_drug.logPow;
logDvow_s = parameters_drug.logDvow-log10(parameters_drug.Y_hh);
parameters_drug.Dvow_s = 10^logDvow_s; 

%% partition coefficients

% Poulin&Thiel
% Berezhkovsky
% 1
part_coeff = 'Berezhkovsky';

[Ptp_na, Ptp_a] = getPartitionCoefficient2(part_coeff, fractions, parameters_drug, [], pH);

Ptp{1} = {Ptp_na, Ptp_a};

%% reorganize pop related parameters

%%% model structure
% 1 - All the organ well stirred
% 2 - PL for every organ
% 3 - Organ PL specified
choose_net = 1;
drug_names = {'MDZ'};

bere_correction = 0;
[ parameters_pbpk2 ] = reorganizeParametersPBPK2( SubjParamOutput , parameters_drug, Ptp , n_samples , drug_names, bere_correction, choose_net, [] );


%% enzyme related parameters

% MPPGL (*11)
MPPGL_mean = 39.79066; % [mg_prot/g_liver] 
MPPGL_CV = 0.269;
MPPGL_sigma = MPPGL_CV*MPPGL_mean;

[ lMPPGL_mean, lMPPGL_sigma ] = getLogNormalValues( 0, MPPGL_mean, MPPGL_sigma );

% CYP microsomal abundance liver (*12)
CYP3A4_mean = 137; %[pmol/mg_microsomal_protein]
CYP3A4_CV = 0.41;
CYP3A4_sd = CYP3A4_mean*CYP3A4_CV;
[ lCYP3A4_mean, lCYP3A4_sd ] = getLogNormalValues( 0, CYP3A4_mean, CYP3A4_sd );

CYP3A5_mean = 103; %[pmol/mg_microsomal_protein]
CYP3A5_CV = 0.65;
CYP3A5_sd = CYP3A5_mean*CYP3A5_CV;
[ lCYP3A5_mean, lCYP3A5_sd ] = getLogNormalValues( 0, CYP3A5_mean, CYP3A5_sd );

corr_l3A4_l3A5 = 0.5228;

%% metabolism related fixed parameters

mw = parameters_drug.mw;
time_corr = 60;

% from (*13)
% Vmax [pmol/min/pmol p450] -> [mg/h/pmol p450]
% Km [uM = u mol/L] -> [mg/L]
Vmax_3a4_1 = 1.96*mw*10^-9*time_corr;
Km_3a4_1 = 2.69*mw*10^-3;
Vmax_3a4_4 = 2.52*mw*10^-9*time_corr;
Km_3a4_4 = 29*mw*10^-3;

Vmax_3a5_1 = 6.7*mw*10^-9*time_corr;
Km_3a5_1 = 10.7*mw*10^-3;
Vmax_3a5_4 = 0.52*mw*10^-9*time_corr;
Km_3a5_4 = 12.1*mw*10^-3;


%% add variability to enzymes parameters


% MPPGL
MPPGL = exp(lMPPGL_sigma*UZ(:,4) + lMPPGL_mean);

% enzymatic abundance

if flag_corr == 1 % correlation
    lambda = sqrt(corr_l3A4_l3A5);
    sigma_eps = sqrt(1 - lambda^2);
    eta = UZ(:,7);
else % no correlation
    lambda = 0;
    sigma_eps = 1;
    eta = 0;
end

eps_3a4 = sigma_eps*UZ(:,5);
eps_3a5 = sigma_eps*UZ(:,6);

lCYP3A4_norm = lambda*eta + eps_3a4;
lCYP3A5_norm = lambda*eta + eps_3a5;

CYP3A4_conc = exp(lCYP3A4_mean + lCYP3A4_sd*lCYP3A4_norm);
CYP3A5_conc = exp(lCYP3A5_mean + lCYP3A5_sd*lCYP3A5_norm);


%% pre allocate structure

P_metabolism(n_samples).MPPGL = [];
P_metabolism(n_samples).CYP3A4_conc = [];
P_metabolism(n_samples).CYP3A5_conc = [];
P_metabolism(n_samples).Vmax_3a4_1 = [];
P_metabolism(n_samples).Km_3a4_1 = [];
P_metabolism(n_samples).Vmax_3a4_4 = [];
P_metabolism(n_samples).Km_3a4_4 = [];
P_metabolism(n_samples).Vmax_3a5_1 = [];
P_metabolism(n_samples).Km_3a5_1 = [];
P_metabolism(n_samples).Vmax_3a5_4 = [];
P_metabolism(n_samples).Km_3a5_4 = [];

%% reorganize parameters & create output cell


for i = 1:n_samples
    
    P_metabolism(i).MPPGL = MPPGL(i);
    P_metabolism(i).CYP3A4_conc = CYP3A4_conc(i);
    P_metabolism(i).CYP3A5_conc = CYP3A5_conc(i);
    
    P_metabolism(i).Vmax_3a4_1 = Vmax_3a4_1;
    P_metabolism(i).Km_3a4_1 = Km_3a4_1;
    P_metabolism(i).Vmax_3a4_4 = Vmax_3a4_4;
    P_metabolism(i).Km_3a4_4 = Km_3a4_4;
    
    P_metabolism(i).Vmax_3a5_1 = Vmax_3a5_1;
    P_metabolism(i).Km_3a5_1 = Km_3a5_1;
    P_metabolism(i).Vmax_3a5_4 = Vmax_3a5_4;
    P_metabolism(i).Km_3a5_4 = Km_3a5_4;
    
end



end


%% references

% (*1)  Gastroplus standard human fasted values
% (*3)  Masoud Jamei, David Turner, Jiansong Yang, Sibylle Neuhoff, Sebastian Polak, Amin Rostami-Hodjegan, and Geoffrey Tucker. Population-Based Mechanistic Prediction of Oral Drug Absorption. The AAPS Journal, Vol. 11, No. 2, June 2009
% (*4)  Lawrence X. Yu, John R. Crison , Gordon L. Amidon (1996). Compartmental transit and dispersion model analysis of small intestinal transit flow in humans. International Journal of Pharmaceutics 140 (1996) 111-118.
% (*5)  N. Parrott, V. Lukacova, G. Fraczkiewicz, and M. B. Bolger. Predicting Pharmacokinetics of Drugs Using Physiologically Based Modeling-Application to Food Effects. The AAPS Journal, Vol. 11, No. 1, March 2009 (2009).
% (*6)  A. Olivares-Morales, Avijit Ghosh, Leon Aarons, and Amin Rostami-Hodjegan. Development of a Novel Simplified PBPK Absorption Model to Explain the Higher Relative Bioavailability of the OROS Formulation of Oxybutynin. The AAPS Journal, Vol. 18, No. 6, November 2016.
% (*8)  FDA drug label, access 14/11/2019 https://www.accessdata.fda.gov/drugsatfda_docs/label/2017/208878Orig1s000lbl.pdf 
% (*9)  PKSim version 7.4
% (*10) Brill et al 2016, doi:10.1002/psp4.12048
% (*11) Simcyp version 16 11/2017
% (*12) Cubitt et al 2011, 
% (*13) Galetin et al 2004 https://doi.org/10.1124/dmd.104.000844 



