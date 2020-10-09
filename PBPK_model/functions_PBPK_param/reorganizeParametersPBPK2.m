function [ parameters_pbpk ] = reorganizeParametersPBPK2( SubjParamOutput , parameters_drug, Ptp , n_samples , drug_name, correction, choose_net, index_organ_PL )

%REORGANIZEPARAMETERSPBPK2 Summary of this function goes here
% union of reorganizeParametersPBPK and berezkowsy correction!

% Nicola Melillo 30/08/2017

%% Output

%%% parameters_pbpk: cell array (n_samples, n_drug)
%   {x}.drug_name: string
%   {x}.V_organs: [L] organs volume, if correction is equal to 1 organs volume
%                 contains blood amount and arterial and venous blood
%                 do not include organs blood content. Otherwhise organs
%                 volume do not contain blood amount that is reparted
%                 between arterial and venous blood.
%   {x}.V_ex_organs: [L] extracellular water volume, same thing as V_organs for
%                    correction.
%   {x}.V_int_organs: [L] intracellular water volume, correction do not
%                     touch it.
%   {x}.Q_organs: [L/min] organs blood flows
%   {x}.V_rotb: [L] rest of the body volume
%   {x}.Q_rotb: [L/min] rest of the body blood flow
%   {x}.Ptp_na: non adipose partition coefficients
%   {x}.Ptp_a: adipose partition coefficients

%% Re-organization of the parameters

sdn = size(drug_name);
n_drug = sdn(1);

parameters_pbpk = cell(n_samples,n_drug);

n_pbpk_eq = length(SubjParamOutput{1}.V_organs);

for i=1:n_drug
    
    Ptp_na_a = Ptp{i};
    Ptp_na = Ptp_na_a{1};
    Ptp_a = Ptp_na_a{2};
    
    for j=1:n_samples
        
        Ptp_na_j = Ptp_na;
        Ptp_a_j = Ptp_a;
        Ptpj = [Ptp_a_j(1);Ptp_na_j(2:end)];
        
        paramPBPKji.drug_name = drug_name{i};
        paramPBPKji.V_iw_organs = SubjParamOutput{j}.V_iw_organs;
        paramPBPKji.Q_organs = SubjParamOutput{j}.Q_organs;
        paramPBPKji.W_organs = SubjParamOutput{j}.W_organs;
        paramPBPKji.Ptp_na = Ptp_na_j;
        paramPBPKji.Ptp_a = Ptp_a_j;
        
        paramPBPKji.V_organs = SubjParamOutput{j}.V_organs;
        paramPBPKji.V_ew_organs = SubjParamOutput{j}.V_ew_organs;
        
        if correction==1
            parameters_drug_ji = parameters_drug{j,i};
            BP = parameters_drug_ji.BP;
            correction_V = SubjParamOutput{j}.V_blood_organs(1:end-2)./Ptpj*BP;
            % correct total V
            paramPBPKji.V_organs(1:end-2) = paramPBPKji.V_organs(1:end-2) + correction_V;
            paramPBPKji.V_organs(end-1:end) = SubjParamOutput{j}.V_blood_organs(end-1:end);
            % correct ew
            paramPBPKji.V_ew_organs = paramPBPKji.V_ew_organs + correction_V;
        end
        
        %%% Volume vector for pbpk system (**1)
        % choose net
        
        if choose_net(i) == 1 
            % 1 - All the organ well stirred
            paramPBPKji.V_pbpk_system = paramPBPKji.V_organs;
            
        elseif choose_net(i) == 2 
            % 2 - PL for every organ
            paramPBPKji.V_pbpk_system = zeros(n_pbpk_eq,1);
            paramPBPKji.V_pbpk_system(1:end-2) =  paramPBPKji.V_ew_organs;
            paramPBPKji.V_pbpk_system(end-1:end) = paramPBPKji.V_organs(end-1:end);
            
        elseif choose_net(i) == 3 
            % 3 - Organ PL specified
            Vpbpk = paramPBPKji.V_organs;
            Vpbpk(index_organ_PL{i}) = paramPBPKji.V_ew_organs(index_organ_PL{i});
            paramPBPKji.V_pbpk_system = Vpbpk;
        end
        
        parameters_pbpk{j,i} = paramPBPKji;
        
    end
    
end



end

%% Comments
% (**1) Why do I create V_pbpk_system? Let's imagine that I have all the
%       organs well-stirred except the liver that is Permeability-Limited. "Distributional PBPK"
%       equation for well-stirred model are equal to these for PL, onique
%       exception for the volume, that in WS model is the total but for PL
%       model is only the extracellular.














