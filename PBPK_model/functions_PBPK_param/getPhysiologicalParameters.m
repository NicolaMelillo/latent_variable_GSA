function [ SubjParamOutput, pH ] = getPhysiologicalParameters( compNames, fractions, par_cell_FM, UZ)

% Nicola Melillo


%% Outputs

%%% SubjParamOutput: cell array, elements are of class ClassSubjParam
%   {x}.W_organs: [kg] organs weight without blood
%   {x}.V_organs: [L] organs volume without blood
%   {x}.Q_organs: [L/min] organs blood flow
%   {x}.W_blood_organs: [kg] blood weight in each organ
%   {x}.V_blood_organs: [L] blood volume in each organ

%%% pH: structure containing pH of various sites
%   .plasma: self explained
%   .iw: intracelluar water of all organs
%   .ew: extracellular water of all organs
%   .bc: blood cells


%% Create Volumes, Blood Volumes and Fluxes

num_comp_av = compNames.num_comp_av;
n_samples = size(UZ,1);

% population parameters (female | male) (*1)
mean_height = [163.3 176.7];  % [cm]
std_height  = [5.85 6.15];    % [cm]
mean_CO     = [6 7] * 60;     % [L/min] -> [L/h]

% Body Mass Index (BMI) ranges for normal weight (*2)
% BMI = weight/height^2, with height in m and weight in kg
BMI = [18.5 24.9]; % [kg/m^2]

% Zd = makedist('Normal', 0, 1);
% Z  = icdf(Zd, U(:,2));

% arterial and venous proportion on total blood (derived from par_cell{i}.WB_organs)
arterial_prop = 1/4;
venous_prop = 3/4;

SubjParamOutput = cell(n_samples, 1);

for i=1:n_samples
    
    % extract the parameters (hp linear scale in function of the height)
    sex_i = double(UZ(i,1) > 0.5); % 0 female | 1 male
    height_i = mean_height(sex_i+1) + UZ(i,2) * std_height(sex_i+1);
    BMI_i = BMI(1) + UZ(i,3) * (BMI(2) - BMI(1));
    
    % derive weight and CO
    weight_i = BMI_i * (height_i/100)^2;
    CO_i = (height_i/mean_height(sex_i+1))^0.75 * mean_CO(sex_i+1); % allometric scaling equation from (*5)
    
    par_cell = par_cell_FM{sex_i+1};
    
    % initialization
    W_organs = zeros(num_comp_av,1);
    Q_organs = zeros(num_comp_av,1);
    V_organs = zeros(num_comp_av,1);
    
    % Organs weight with blood
    W_organs_with_blood = weight_i * par_cell.W_f;
    W_blood_organs = W_organs_with_blood(end) * par_cell.WB_f;
    V_blood_organs = W_blood_organs; % hp: density=1 [kg/L]
    
    % organs volumes excluded arterial and venous blood
    % first remove blood content from organ weight, pay attention on what
    % weight fractions are you using! Willmann contain blood content, other
    % may not contain it (e.g. from ICRP).
    W_organs(1:end-2) = W_organs_with_blood(1:end-1) - W_blood_organs(1:end-2);
    V_organs(1:end-2) = W_organs(1:end-2)./fractions.density; % [kg]*([kg/L])^-1 -> [L] (**1)
    Q_organs(1:end-2) = CO_i * par_cell.Q_f;
    
    % Total blood volume parted between arterial and venous blood (**2)
    W_organs(end-1) = W_organs_with_blood(end)*arterial_prop;
    W_organs(end) = W_organs_with_blood(end)*venous_prop;
    V_organs(end-1) = W_organs(end-1); % hp: density=1 [kg/L]
    V_organs(end) = W_organs(end); % hp: density=1 [kg/L]

    % Intracellular and extracellular volume
    V_iw_organs = V_organs(1:end-2).*fractions.Fiw_RR;
    V_ew_organs = V_organs(1:end-2).*fractions.Few_RR;
    
    % output
    SubjParam_i.sex = sex_i;
    SubjParam_i.height = height_i;
    SubjParam_i.weight = weight_i;
    SubjParam_i.CO = CO_i;
    SubjParam_i.W_organs = W_organs; % [kg]
    SubjParam_i.V_organs = V_organs; % [L]
    SubjParam_i.Q_organs = Q_organs; % [L/h]
    SubjParam_i.W_blood_organs = W_blood_organs; % [kg]
    SubjParam_i.V_blood_organs = V_blood_organs; % [L]
    SubjParam_i.V_iw_organs = V_iw_organs; % [L]
    SubjParam_i.V_ew_organs = V_ew_organs; % [L]
    
    SubjParamOutput{i} = SubjParam_i;
    
end

%% pH in various parts of the body

pH.plasma = 7.4;
pH.ew = 7.4; % (*3) in general is the extracellular water pH of all organs
pH.iw = 7; % (*3) intracellular water pH
pH.bc = 7.22; % (*4) blood cells pH


end


%% References

% (*1) Cacciari 2006, Journal of Endocrinological Investigation, https://doi.org/10.1007/BF03344156
% (*2) WHO BMI ranges for "Normal weight" class https://www.euro.who.int/en/health-topics/disease-prevention/nutrition/a-healthy-lifestyle/body-mass-index-bmi
% (*3) TRUDY RODGERS, MALCOLM ROWLAND 2005 Physiologically Based Pharmacokinetic Modelling 2:Predicting the Tissue Distribution of Acids,Very Weak Bases, Neutrals and Zwitterions
% (*4) Trudy Rodgers and Malcolm Rowland Mechanistic Approaches to Volume of Distribution Predictions: Understanding the Processes 2007
% (*5) Willmann 2007 JPKPD https://doi.org/10.1007/s10928-007-9053-5

%% Comments

% (**1) fractions.density is the specific gravity of each tissue, ergo
%       density_tissue/density_water. Given that density_water is equal to 1
%       [kg/L] so density_tissue = fractions.density*density_water that,
%       except the measure unit, it is equal to fractions.density.

% (**2) Blood content is distributed in: organs, large arteries and large
%       veins. In W_organs (and V_organs) I consider for each organ the
%       volume WITHOUT blood. Arterial and venous volumes (and weight) are
%       respectively 1/4 and 3/4 of TOTAL blood volume













