function [Y] = PBPK_uqlab(U, flg, flg_gen)

%% load parameters

n_samples = size(U,1);

[P_pbpk, P_drug, P_metabolism] = load_physio_param(U,flg);

cyp3a5_gen = ones(n_samples,1) * flg_gen;


%% execute the model

dose = 5; % [mg]
tspan = [0 24*7];

n_eq_pbpk = 17;

X0 = zeros(n_eq_pbpk,1);
X0(end) = dose;

Y = zeros(n_samples, 1);

parfor i = 1:n_samples
    
    dX_c = @(t,X) midazolam_pbpk(t, X, P_pbpk{i}, P_drug, P_metabolism(i), cyp3a5_gen(i));
    [t_c, X_c] = ode15s(dX_c, tspan, X0);
    
    V_ven = P_pbpk{i}.V_pbpk_system(end);
    Y(i) = trapz(t_c, X_c(:,end))./V_ven/P_drug.BP; % plasma
    
end

end

