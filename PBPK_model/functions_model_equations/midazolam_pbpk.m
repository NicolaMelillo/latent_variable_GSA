function [dX] = midazolam_pbpk(t, X, P_pbpk, P_drug, P_metabolism, cyp3a5_gen)

%% info

% 1 'adipose';
% 2 'bone';
% 3 'brain';
% 4 'heart';
% 5 'muscle';
% 6 'skin';
% 7 'spleen';
% 8 'kidney';
% 9 'gonads';
% 10 'lung';
% 11 'stomach';
% 12 'S int';
% 13 'L int';
% 14 'liver';
% 15 'pancreas';
% 16 'arterial';
% 17 'venous';

%% parameters

% ------------- physiological parameters
V_liv = P_pbpk.V_organs(14); % [L]
W_liv = P_pbpk.W_organs(14)*1000; % [kg] -> [g]

MPPGL = P_metabolism.MPPGL;
CYP3A4_conc = P_metabolism.CYP3A4_conc;
CYP3A5_conc = P_metabolism.CYP3A5_conc;

% ------------- drug related parameters

Vmax_3a4_1 = P_metabolism.Vmax_3a4_1*MPPGL*W_liv*CYP3A4_conc;
Km_3a4_1 = P_metabolism.Km_3a4_1*V_liv;
Vmax_3a4_4 = P_metabolism.Vmax_3a4_4*MPPGL*W_liv*CYP3A4_conc;
Km_3a4_4 = P_metabolism.Km_3a4_4*V_liv;

Vmax_3a5_1 = P_metabolism.Vmax_3a5_1*MPPGL*W_liv*CYP3A5_conc;
Km_3a5_1 = P_metabolism.Km_3a5_1*V_liv;
Vmax_3a5_4 = P_metabolism.Vmax_3a5_4*MPPGL*W_liv*CYP3A5_conc;
Km_3a5_4 = P_metabolism.Km_3a5_4*V_liv;

%% transformed parameters

X_liv = X(14) * P_drug.fut;

MET_3a4 = Vmax_3a4_1/(Km_3a4_1 + X_liv) + Vmax_3a4_4/(Km_3a4_4 + X_liv);
MET_3a5 = Vmax_3a5_1/(Km_3a5_1 + X_liv) + Vmax_3a5_4/(Km_3a5_4 + X_liv);
MET_3a5 = MET_3a5 * cyp3a5_gen;

MET_h  = (MET_3a4 + MET_3a5) * X_liv; % [mg/min]

%% equation system

% distribution pbpk
[ dX_pbpk ] = pbpk_transport_02( t, X, P_pbpk, P_drug.BP );

% add metabolism within the liver
dX_pbpk(14) = dX_pbpk(14) - MET_h;

dX = dX_pbpk;

end

