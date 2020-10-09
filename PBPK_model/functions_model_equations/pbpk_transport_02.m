function [ dA ] = pbpk_transport_02( t, A, parameters_pbpk, BP )

% Represent only transport functions in pbpk, no CL, no network, no input
% from bolus or any ADAM model.
% Only raw structure.

% C in equations MUST be a concentration, A in argument is an amount, so I
% have to impose C = A/Volumes

% for well stirred model delete PL_par in inputs and select well stirrend
% in liver part.

% Compartments article order / Compartments my order
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

%% Parameters

V = parameters_pbpk.V_pbpk_system;
Q = parameters_pbpk.Q_organs;
Ptp_na = parameters_pbpk.Ptp_na;
Ptp_a = parameters_pbpk.Ptp_a;

% system charateristics
n_eq_pbpk = length(V);
dC = zeros(n_eq_pbpk,1);
i = 1;

%% Equations

% set concentration
C = A./V;

C_art = C(end-1);
C_ven = C(end);

% 1-Adipose
dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_a(i)/BP) );
venous_eq = Q(i)*C(i)/(Ptp_a(i)/BP);
arterial_eq = -Q(i)*C_art;

% from comp 2 to comp 6 my order (bone,brain,heart,muscle,skin)
for i=2:6
   dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_na(i)/BP) );
   venous_eq = venous_eq+Q(i)*C(i)/(Ptp_na(i)/BP);
   arterial_eq = arterial_eq-Q(i)*C_art;
end

% 7-Spleen
i=i+1;
dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_na(i)/BP) );
arterial_eq = arterial_eq-Q(i)*C_art;

% 8-Kidney
i = i+1;
dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_na(i)/BP) );
venous_eq = venous_eq+Q(i)*C(i)/(Ptp_na(i)/BP);
arterial_eq = arterial_eq-Q(i)*C_art;

% 9-Gonads
i = i+1;
dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_na(i)/BP) );
venous_eq = venous_eq+Q(i)*C(i)/(Ptp_na(i)/BP);
arterial_eq = arterial_eq-Q(i)*C_art;

% 10-Lung
i = i+1;
dC(i) = Q(i)/V(i)*( C_ven - C(i)/(Ptp_na(i)/BP));
venous_eq = venous_eq - Q(i)*C_ven;
arterial_eq = arterial_eq+Q(i)*C(i)/(Ptp_na(i)/BP);

% 11-Stomach
i = i+1;
dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_na(i)/BP) );
arterial_eq = arterial_eq-Q(i)*C_art;

% 12-S_int
i = i+1;
dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_na(i)/BP) );
arterial_eq = arterial_eq-Q(i)*C_art;

% 13-L_int
i = i+1;
dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_na(i)/BP) );
arterial_eq = arterial_eq-Q(i)*C_art;

% 14-Liver

% portal input for liver
portal_index = [7, 11, 12, 13, 15];
lp = length(portal_index);
portal_eq = 0;
Q_port_ven = 0;
for j=1:lp
    jj = portal_index(j);
    portal_eq = portal_eq + Q(jj)*( C(jj)/(Ptp_na(jj)/BP) );
    Q_port_ven = Q_port_ven + Q(jj);
end

i = i+1;
Q_art_liver = Q(i);
Q_ven_liver = Q_port_ven + Q_art_liver;
dC(i) = ( ( Q_art_liver )*C_art + portal_eq - Q_ven_liver*C(i)/(Ptp_na(i)/BP) )/V(i); %well stirred
venous_eq = venous_eq + Q_ven_liver*C(i)/(Ptp_na(i)/BP);
arterial_eq = arterial_eq - ( Q_art_liver )*C_art;

% 15-Pancreas
i = i+1;
dC(i) = Q(i)/V(i)*( C_art - C(i)/(Ptp_na(i)/BP) );
arterial_eq = arterial_eq-Q(i)*C_art;

% 16-Arterial Blood
i = i+1;
dC(i) = arterial_eq/V(i);

% 17-Venous Blood
i = i+1;
dC(i) = venous_eq/V(i);

dA = dC.*V; % output in amount


end

