function [ Ptp_na, Ptp_a ] = getPartitionCoefficient2( type, fractions, parameters_drug, SubjOutput, pH )
%GETPARTITIONCOEFFICIENT2 To obtain partition coefficient from drug
%information and subjects info.

% Nicola Melillo 28/08/2017

% common drug parameters
Pow = parameters_drug.Pow;
fup = parameters_drug.fup;
Dvow_s = parameters_drug.Dvow_s;
Dvow = 10^parameters_drug.logDvow;
type_drug = parameters_drug.type;


switch type
    
    case 'Poulin&Thiel' % (*1)
%%        
        % drug parameters
        fut = parameters_drug.fut;
        
        % input volume parameters Poulin
        Vw = fractions.Vw_poulin;
        Vnl = fractions.Vnl_poulin;
        Vph = fractions.Vph_poulin;
        
        % rearranging parameters
        Vw_t = Vw(1:end-1);
        Vw_p = Vw(end);
        Vnl_t = Vnl(1:end-1);
        Vnl_p = Vnl(end);
        Vph_t = Vph(1:end-1);
        Vph_p = Vph(end);
        
        % partition coefficients
        Ptp_na = (  Dvow*(Vnl_t + 0.3*Vph_t) + (Vw_t + 0.7*Vph_t) )./( Dvow*(Vnl_p + 0.3*Vph_p) + (Vw_p + 0.7*Vph_p) )*fup/fut; % non adipose partition coefficient tissue:plasma
        Ptp_a = (  Dvow_s*(Vnl_t + 0.3*Vph_t) + (Vw_t + 0.7*Vph_t) )./( Dvow_s*(Vnl_p + 0.3*Vph_p) + (Vw_p + 0.7*Vph_p) )*fup; % adipose partition coefficient
        
    case 'Berezhkovsky' % (*2)
%%        
        % drug parameters
        fut = parameters_drug.fut;
        
        % input volume parameters Poulin
        Vw = fractions.Vw_poulin;
        Vnl = fractions.Vnl_poulin;
        Vph = fractions.Vph_poulin;
        
        % rearranging parameters
        Vw_t = Vw(1:end-1);
        Vw_p = Vw(end);
        Vnl_t = Vnl(1:end-1);
        Vnl_p = Vnl(end);
        Vph_t = Vph(1:end-1);
        Vph_p = Vph(end);
        
        
        % partition coefficients
        Ptp_a = (  Dvow*(Vnl_t + 0.3*Vph_t) + (Vw_t/fut + 0.7*Vph_t) )./( Dvow*(Vnl_p + 0.3*Vph_p) + (Vw_p/fup + 0.7*Vph_p) ); % adipose partition coefficient tissue:plasma
        Ptp_na = Ptp_a;
        
    case '1'
%%  
        AR_tp = fractions.AR_tp_RR; % prima c'era LP_tp
        l = length(AR_tp);
        Ptp_na = 1*ones(l,1);
        Ptp_a = 1*ones(l,1);
        
end

end


%% references partition coefficients
% (*1) PATRICK POULIN, FRANK-PETER THEIL 2001 Prediction of Pharmacokinetics Prior to In Vivo Studies.1. Mechanism-Based Prediction of Volume of Distribution
% (*2) LEONID M. BEREZHKOVSKIY 2004 Volume of Distribution at Steady State for a Linear Pharmacokinetic System with Peripheral Elimination, Journal of Pharmaceutical Sciences


