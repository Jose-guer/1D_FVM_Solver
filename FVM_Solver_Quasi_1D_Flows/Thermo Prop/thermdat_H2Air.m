function [gamma, enthalpy, R, cp_mix] = thermdat_H2Air(Trange,Y,mol_units)
%Inputs :
% Trange [N x 1] - temperature [K]
% Y [N,4] - [YH2O,YH2,YN2,YO2] mass fractions
% 
% Y = 
% [YH2O(T1), YH2(T1), YN2(T1), YO2(T1);
% [YH2O(T2), YH2(T2), YN2(T2), YO2(T2);
% .
% .
% .
% [YH2O(TN), YH2(TN), YN2(TN), YO2(TN);

%Outputs :
% gamma - specific heats ratio of mixture
% enthalpy of each component at T [J/kg]


% mixture file must have [x], [poly], [mw], and [visco] defined such that:
%x = [molefrac1 molefrac2...]
%poly = [polynomial(poly#,Trange,component)]
%mw = [component1 component2 ...] %kg/mol
%polynomial fit is A*T^-2+B*T^-1+C+D*T+E*T^2+F*T^3+G*T^4
%visco = [c To etao]

% ***** NASA9 Polynomial Data Documentation *****
% McBride, Bonnie J. NASA Glenn coefficients for calculating thermodynamic
% properties of individual species. National Aeronautics and Space
% Administration, John H. Glenn Research Center at Lewis Field, 2002.

% ***** Data Format *****
%Note that H0 and S0 have integration constants b1 and b2
%H = H(298.15K) + [H(T) - H(298.15K)]

% Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4

% H/RT = -a1*T^-2 + a2*T^-1*ln(T) + a3 + a4*T/2 + a5*T^2/3 +
%                                               a6*T^3/4 + a7*T^4/5 + b1/T

% S/R = -a1*T^-2/2 - a2*T^-1 + a3*ln(T) + a4*T + a5*T^2/2 +
%                                                  a6*T^3/6 + a7*T^4/4 + b2

%Data source
% https://shepherd.caltech.edu/EDL/PublicResources/sdt/SDToolbox/cti/NASA9/nasa9.dat
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

%% NASA-9 Polynomial Data

R = 8.3145; %J/mol-K

species = {'H2O','H2','N2','O2'};
N_species = 4;
L = length(Trange);

MW_arr = zeros(1,N_species);
cp = zeros(L,N_species);
enthalpy = zeros(L,N_species);

for j = 1:N_species

    comp = species{j};

    switch comp

        case 'H2O'

            poly = [-3.947960830e+04 5.755731020e+02 9.317826530e-01 7.222712860e-03 -7.342557370e-06 4.955043490e-09 -1.336933246e-12  -3.303974310e+04 1.724205775e+01;
                1.034972096e+06 -2.412698562e+03 4.646110780e+00 2.291998307e-03 -6.836830480e-07 9.426468930e-11 -4.822380530e-15  -1.384286509e+04 -7.978148510e+00 ];

            MW = 18.0153;

        case 'H2'

            poly = [4.078322810e+04 -8.009185450e+02 8.214701670e+00 -1.269714360e-02 1.753604930e-05 -1.202860160e-08 3.368093160e-12  2.682484380e+03 -3.043788660e+01;
                5.608123380e+05 -8.371491340e+02 2.975363040e+00 1.252249930e-03 -3.740718420e-07 5.936628250e-11 -3.606995730e-15   5.339815850e+03 -2.202764050e+00;
                4.966716130e+08 -3.147448120e+05 7.983887500e+01 -8.414504190e-03 4.753060440e-07 -1.371809730e-11 1.605374600e-16 2.488354660e+06 -6.695524190e+02];

            MW = 2.016;

        case 'O2'

            poly = [ -3.425563420e+04 4.847000970e+02 1.119010961e+00 4.293889240e-03 -6.836300520e-07 -2.023372700e-09 1.039040018e-12 -3.391454870e+03 1.849699470e+01;
                -1.037939022e+06 2.344830282e+03 1.819732036e+00 1.267847582e-03 -2.188067988e-07 2.053719572e-11 -8.193467050e-16 -1.689010929e+04 1.738716506e+01;
                4.975294300e+08 -2.866106874e+05 6.690352250e+01 -6.169959020e-03 3.016396027e-07 -7.421416600e-12 7.278175770e-17 2.293554027e+06 -5.530621610e+02];

            MW = 31.999;

        case 'N2'

            poly = [2.210371497E+04,-3.818461820E+02,6.082738360E+00,-8.530914410E-03,1.384646189E-05,-9.625793620E-09,2.519705809E-12,7.108460860E+02,-1.076003744E+01;
                5.877124060E+05,-2.239249073E+03,6.066949220E+00,-6.139685500E-04,1.491806679E-07,-1.923105485E-11,1.061954386E-15,1.283210415E+04,-1.586640027E+01;
                8.310139160E+08,-6.420733540E+05,2.020264635E+02,-3.065092046E-02,2.486903333E-06,-9.705954110E-11,1.437538881E-15,4.938707040E+06,-1.672099740E+03];

            MW = 28.0134;
    end


    MW_arr(j) = MW;

    h = @(coeff,T,R) -coeff(1,1)*R*T^-1 +...
        coeff(1,2)*R*log(T) +...
        coeff(1,3)*R*T +...
        coeff(1,4)*R*(T^2)/2 +...
        coeff(1,5)*R*(T^3)/3 +...
        coeff(1,6)*R*(T^4)/4 +...
        coeff(1,7)*R*(T^5)/5 +...
        coeff(1,8)*R;

    for i = 1:L
        T = Trange(i);

        if T < 1000
            cp(i,j)= R*( poly(1,1)*T^-2 + poly(1,2)*T^-1 + poly(1,3) + poly(1,4)*T + poly(1,5)*T^2 + poly(1,6)*T^3 + poly(1,7)*T^4 ); %J/mol-K
            enthalpy(i,j) = h(poly(1,:),T,R); %J/mol

        elseif T < 6000
            cp(i,j)= R*( poly(2,1)*T^-2 + poly(2,2)*T^-1 + poly(2,3) + poly(2,4)*T + poly(2,5)*T^2 + poly(2,6)*T^3 + poly(2,7)*T^4 ); %J/mol-K
            enthalpy(i,j) = h(poly(2,:),T,R); %J/mol

        else
            cp(i,j)= R*( poly(3,1)*T^-2 + poly(3,2)*T^-1 + poly(3,3) + poly(3,4)*T + poly(3,5)*T^2 + poly(3,6)*T^3 + poly(3,7)*T^4 ); %J/mol-K
            enthalpy(i,j) = h(poly(3,:),T,R); %J/mol
        end

    end

end


% i = 1, [h(H2O,T1), h(H2,T1), h(N2,T1), h(O2,T1);
% i = 2, [h(H2O,T2), h(H2,T2), h(N2,T2), h(O2,T2);
% .
% .
% .
% i = N, [h(H2O,TN), h(H2,TN), h(N2,TN), h(O2,TN);

%% Convert to units

if mol_units == 1

    MW_mix = (sum(Y./MW_arr,2)).^(-1); %kg/kmol

    X = Y .* MW_mix ./ MW_arr;

    enthalpy = enthalpy * 1000;  %J/kmol-K
    cp = cp * 1000;              %J/kmol-K
    R = R * 1000;                %J/kmol-K

    cp_mix = sum(cp.*X,2);       %J/kmol-K
    gamma  = cp_mix./(cp_mix-R); %J/kmol-K


else
    
    MW_arr = MW_arr./1000;             %kg/mol
    enthalpy = enthalpy ./ MW_arr;     %J/kg
    cp = cp./MW_arr;                   %J/kg-K

    MW_mix = (sum(Y./MW_arr,2)).^(-1); %kg/kmol
    cp_mix = sum(cp.*Y,2);             %J/kg-K
    R = R./MW_mix;                     %J/kg-K
    gamma  = cp_mix./(cp_mix-R);       %J/kg-K

    % enthalpy_mix = sum(enthalpy.*Y,2); %For the Mixture , J/kg

end

end


