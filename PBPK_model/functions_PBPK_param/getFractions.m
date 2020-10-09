function [compNames, fractions, par_cell_FM ] = getFractions()

% Nicola Melillo, UoM 24/09/2020

%% Outputs

% All output, except scalars and the first, are classes

%%% compNames (structure);
%   .compartments_my_order: all organs & tissues names, the last is blood
%   .compartments_my_order_wo_blood: all organs & tissues names except blood
%   .compartments_my_order_av: all organs & tissues names, the last two are arterial and venous blood
%   .mapCompNames: hashmap that links organs name to the index in compartments_my_order
%   .num_comp_av: number of organs including arterial and venous blood

%%% fractions (ClassFraction);
%   .Vw_poulin: with blood (plasma) Poulin
%   .Vnl_poulin: with blood (plasma) Poulin
%   .Vph_poulin: with blood (plasma) Poulin
%   .Fnl_RR: with blood (plasma) Rogertz & Rowland
%   .Fnp_RR: with blood (plasma) Rogertz & Rowland
%   .Fw_RR: with blood (plasma) Rogertz & Rowland
%   .Fiw_RR: without blood (plasma) Rogertz & Rowland
%   .Few_RR: without blood (plasma) Rogertz & Rowland
%   .AR_tp_RR: without blood (plasma) Rogertz & Rowland
%   .LP_tp_RR: without blood (plasma) Rogertz & Rowland
%   .AP_t: without blood (plasma) Rogertz & Rowland
%   .AP_tbc: without blood (plasma) Rogertz & Rowland
%   .Fnp_bc_RR: blood cell neutral phospholipids
%   .Fnl_bc_RR: blood cell neutral lipids
%   .tw_bc_RR: blood cell total water
%   .density: [kg/ml]

%%% par_cell_FM: cell containing {par_frac_female, par_frac_male}
%   {x}.W_f: organs weight fractions
%   {x}.Q_f: fraction of cardiac output directed in each organ
%   {x}.WB_f: fraction of blood weight (on total blood weight) present in each organ (also in arterial and venous blood!)



%% Compartments Names

%%% pbpk compartments names
% Poulin & Theil and RR
compartments_article_order = {'adipose';'bone';'brain';'stomach';'S int';'L int';'heart';'kidney';'liver';'lung';'muscle';'pancreas';'skin';'spleen';'gonads';'blood'};
compartments_article_order_bf = {'adipose';'bone';'brain';'stomach';'S int';'L int';'heart';'kidney';'liver';'lung';'muscle';'pancreas';'skin';'spleen';'gonads'};
% ICRP order
parameters_ICRP_bvolume = {'adipose','brain','stomach','S int','L int','heart','kidney','liver','lung','muscle','pancreas','bone','skin','spleen','gonads','venous','arterial'}';
% my order
compartments_my_order = {'adipose';'bone';'brain';'heart';'muscle';'skin';'spleen';'kidney';'gonads';'lung';'stomach';'S int';'L int';'liver';'pancreas';'blood'};
compartments_my_order_wo_blood = {'adipose';'bone';'brain';'heart';'muscle';'skin';'spleen';'kidney';'gonads';'lung';'stomach';'S int';'L int';'liver';'pancreas'};
compartments_my_order_av = {'adipose';'bone';'brain';'heart';'muscle';'skin';'spleen';'kidney';'gonads';'lung';'stomach';'S int';'L int';'liver';'pancreas';'arterial';'venous'};

sc = size(compartments_my_order_av);

% Index hashmap
s_param_tot = size(compartments_my_order);
n_compartments = s_param_tot(1);
mapCompNames = containers.Map(compartments_my_order,1:n_compartments);

% output
compNames.compartments_my_order = compartments_my_order;
compNames.compartments_my_order_wo_blood = compartments_my_order_wo_blood;
compNames.compartments_my_order_av = compartments_my_order_av;
compNames.mapCompNames = mapCompNames;
compNames.num_comp_av = sc(1);

%% Physical composition organs - fractions

% All fraction information for gonads where taken from (*7)
% If information of some organs (e.g. plasma for Fw in RR) in RR was not
% provided I integrated with information found in Poulin and viceversa.

% Fractions: article order Poulin (*2)
% Pancreas from RR (*3) (*4)
Vw_a = [0.18 0.439 0.77 0.718 0.718 0.718  0.758 0.783 0.751 0.811 0.76 0.641 0.718 0.788 0.8 0.945]'; % water fraction (*1)
Vnl_a = [0.79 0.074 0.051 0.0487 0.0487 0.0487 0.0115 0.0207 0.0348 0.003 0.0238 0.0403 0.0284 0.0201 4.8e-3 0.0035]'; % neutral lipids volume (*1)
Vph_a = [0.002 0.0011 0.0565 0.0163 0.0163 0.0163 0.0166 0.0162 0.0252 0.009 0.0072 0.0090 0.0111 0.0198 0.01 0.00225]'; % phospholipids volume (*1)

% Fractions: article order Rogertz, Rowland (*3) (*4)
Fnl_a = [0.0016 0.0174 0.0391 0.0375 0.0375 0.0375 0.0135 0.0121 0.0135 0.0215 0.01 0.0403 0.0603 0.0071 4.8e-3 0.0023]';
Fnp_a = [0.853 0.0016 0.0015 0.0124 0.0124 0.0124 0.0106 0.0240 0.0238 0.0123 0.0072 0.0090 0.0044 0.0107 0.01 0.0013]';
Fw_a = [0.144 0.417 0.753 0.738 0.738 0.738 0.568 0.672 0.642 0.574 0.726 0.658 0.0641 0.562 0.8 0.945]';
Few_a = [0.135 0.1 0.162 0.282 0.282 0.282 0.320 0.273 0.161 0.336 0.118 0.120 0.382 0.207 0.06]';
Fiw_a = [0.017 0.346 0.620 0.475 0.475 0.475 0.456 0.483 0.573 0.446 0.630 0.664 0.291 0.579 0.88]';
AR_tp_a = [0.049 0.1 0.048 0.158 0.158 0.158 0.157 0.130 0.086 0.212 0.064 0.060 0.277 0.097 0.05]'; % albumin tissue to plasma ratio
LP_tp_a = [0.068 0.05 0.041 0.141 0.141 0.141 0.160 0.137 0.161 0.168 0.059 0.060 0.096 0.207 0.04]'; % lipoprotein tissue to plasma ratio
AP_t_a = [0.4 0.67 0.4 2.41 2.41 2.41 2.25 5.03 4.56 3.91 1.53 1.67 1.32 3.18 2.45]'; % acidic phospholipids amount tissue [mg/g]
AP_bc = 0.5;% acidic Phospholipids blood cell

% (*5)
% HP: gonads specific gravity was fixed to 1
specific_gravity_a = [ 0.916 1.4303 1.0365 1.046 1.046 1.046 1.030  1.050 1.08 1.0515 1.041 1.045 1.1754 1.054 1]'; % adimensional

%%% Sort vectors and create output!

% my order Poulin
fractions.Vw_poulin = changeorder_063(compartments_article_order,compartments_my_order,Vw_a); % with blood (plasma)
fractions.Vnl_poulin = changeorder_063(compartments_article_order,compartments_my_order,Vnl_a); % with blood (plasma)
fractions.Vph_poulin = changeorder_063(compartments_article_order,compartments_my_order,Vph_a); % with blood (plasma)

% my order Rogertz, Rowland
fractions.Fnl_RR = changeorder_063(compartments_article_order,compartments_my_order,Fnl_a); % with blood (plasma)
fractions.Fnp_RR = changeorder_063(compartments_article_order,compartments_my_order,Fnp_a); % with blood (plasma)
fractions.Fw_RR = changeorder_063(compartments_article_order,compartments_my_order,Fw_a); % with blood (plasma)
fractions.Fiw_RR = changeorder_063(compartments_article_order_bf,compartments_my_order_wo_blood,Fiw_a); % without blood (plasma)
fractions.Few_RR = changeorder_063(compartments_article_order_bf,compartments_my_order_wo_blood,Few_a); % without blood (plasma)
fractions.AR_tp_RR = changeorder_063(compartments_article_order_bf,compartments_my_order_wo_blood,AR_tp_a); % without blood (plasma)
fractions.LP_tp_RR = changeorder_063(compartments_article_order_bf,compartments_my_order_wo_blood,LP_tp_a); % without blood (plasma)
fractions.AP_t = changeorder_063(compartments_article_order_bf,compartments_my_order_wo_blood,AP_t_a); % without blood (plasma)
fractions.AP_tbc = fractions.AP_t/AP_bc;  % without blood (plasma)
fractions.Fnp_bc_RR = 0.0033; %(*6) blood cell neutral phospholipids
fractions.Fnl_bc_RR = 0.0012; %(*6) blood cell neutral lipids
fractions.tw_bc_RR = 0.63; %(*6) blood cell total water

% my order specific gravity
specific_gravity = changeorder_063(compartments_article_order_bf,compartments_my_order_wo_blood,specific_gravity_a);  % without blood (plasma)
density_water = 1; % [kg/L]
fractions.density = specific_gravity*density_water; % [kg/L]

%% Organs weight and blood flow fractions

% HP: they are constant within the population

% Standard subject ICRP (*8)
Wmale_mean_ICRP = 73; % [kg]
Wfemale_mean_ICRP = 60; % [kg]

% Blood volumes fraction (*8)
% to get the fraction I divided the blood volume for the standard subject
% from the same source
blood_perc_male = 5.600/Wmale_mean_ICRP*100;
blood_perc_female = 4.100/Wfemale_mean_ICRP*100;

% Volumes fractions include blood content (*9) not to order
Wmale_f = [20.4 16.2 2.1 0.57 44.3 5.2 0.33 0.6 0.06 1.8 0.23 1 0.56 3.2 0.26 blood_perc_male]'/100;
Wfemale_f = [32.2 15.2 2.3 0.55 33.8 4.5 0.37 0.67 0.02 1.7 0.27 1.2 0.69 3.2 0.28 blood_perc_female]'/100;

% Blood percentage in organs, relative to total blood in body (*8) to order
Wfrac_blood_male = [5.0 1.2 1.0 3.8 2.2 1.0 2.0 10 10.5 14 0.6 7.0 3.0 1.4 0.04 6.0 18]'/100; % fractions
Wfrac_blood_female = [8.5 1.2 1.0 3.8 2.2 1.0 2.0 10 10.5 10.5 0.6 7.0 3.0 1.4 0.02 6 18]'/100; % fractions
Wfrac_blood_male = changeorder_063(parameters_ICRP_bvolume, compartments_my_order_av, Wfrac_blood_male);
Wfrac_blood_female = changeorder_063(parameters_ICRP_bvolume, compartments_my_order_av, Wfrac_blood_female);

% Blood flow Willmann (*9) not to order
Qperc_male = [5.3 5.3 12.8 4.3 18.1 5.3 3.2 21.7 0.05 100 1.1 10.6 4.3 6.9 1.1]'/100;
Qperc_female = [9 5 13 5 12 5 3 20 0.02 100 1 12 5 7 1]'/100;


%%%
% M
par_frac_male.W_f = Wmale_f;
par_frac_male.Q_f = Qperc_male;
par_frac_male.WB_f = Wfrac_blood_male;
% F
par_frac_female.W_f = Wfemale_f;
par_frac_female.Q_f = Qperc_female;
par_frac_female.WB_f = Wfrac_blood_female;

% output
par_cell_FM{1} = par_frac_female;
par_cell_FM{2} = par_frac_male;


end


%% References

% (*1) PATRICK POULIN, FRANK-PETER THEIL 2001 Prediction of Pharmacokinetics Prior to In Vivo Studies.1. Mechanism-Based Prediction of Volume of Distribution
% (*2) PATRICK POULIN, FRANK-PETER THEIL 2002 Prediction of Pharmacokinetics prior to In Vivo Studies.II. Generic Physiologically Based Pharmacokinetic Models of Drug Disposition
% (*3) TRUDY RODGERS, MALCOLM ROWLAND 2005 Physiologically Based Pharmacokinetic Modelling 2:Predicting the Tissue Distribution of Acids,Very Weak Bases, Neutrals and Zwitterions
% (*4) TRUDY RODGERS, DAVID LEAHY, MALCOLM ROWLAND 2004 Physiologically Based Pharmacokinetic Modeling 1:Predicting the Tissue Distribution of Moderate-to-Strong Bases
% (*5) R.P. Brown, M.D. Delp, L. Lindstedt et al. 1997 PHYSIOLOGICAL PARAMETER VALUES FOR PHYSIOLOGICALLY BASED PHARMACOKINETIC MODELS, Toxicology and Industrial Health
% (*6) Trudy Rodgers1,2 and Malcolm Rowland Mechanistic Approaches to Volume of Distribution Predictions: Understanding the Processes 2007
% (*7) Open System Pharmacology Suite 7.1 release, downloaded in June 2017
% (*8) ICRP 89
%       - adipose is "separable adipose tissue, excluding yellow marrow"
%       - Stomach and intestine weights are Wall weights, not content
%       - Large intestine is the sum of right, left colon plus rectosigmoid
%       - Hearth is "tissue only"
%       - muscle is skeletal
%       - respiratory tissue is lung - tissue only
%       - skeleton is total skeleton
%       - gonads for male are testes, epididymes and prostate
%       - gonads for female are ovaries, fallopian tubes and uterus
% (*9) Willmann, Hohn, Edginton, Sevestre, Solodenko, Weiss, Lippert and Schmitt, "Development of a Physiology-Based Whole-Body Population Model for Assessing the Influence of Individual Variability on the Pharmacokinetics of Drugs", Journal of pharmacokinetics and pharmacodynamics, vol 34, no 3, June 2007
% (*10) Lynn R. Williams "Reference values for total blood volume and cardiac output in humans", 1994



