% This script uses Cantera to compare the thermodynamic properties of
% H2-H2O-N2-O2 mixtures with the outputs of the function thermdat_H2Air.m
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

close all;
clear; clc;

dx = 0.05;
x = 0:dx:1;

phi = 0.6;
etac = x';

[Y,dY,X,dX] = species_var(etac,phi,dx);

% figure; plot(x,Y)

Trange = linspace(1000,2200,length(x));
Trange = Trange';

mol_units = 1;
[gamma, enthalpy,~,cp] = thermdat_H2Air(Trange,Y,mol_units);

gas = Solution('mevel2017.yaml');

h_arr = 0*Trange;
cp_arr = 0*Trange;
gam_arr = 0*Trange;


for j = 1:length(Trange)
    T = Trange(j);

    %Individual properties
    % q = 'H2O:1';

    %Mixture properties
    q = ['H2O:', num2str(Y(j,1)), ' H2:',num2str(Y(j,2)),  ' N2:',num2str(Y(j,3)),   ' O2:',num2str(Y(j,4))    ];
    set(gas,'T',T,'P',oneatm,'MassFractions',q)
    
    %Mass Units
    % cp_arr(j) = cp_mass(gas);
    % h_arr(j) = enthalpy_mass(gas);
    % gam_arr(j) = cp_mass(gas)/cv_mass(gas);

    %Mole Units
    cp_arr(j) = cp_mole(gas);
    h_arr(j) = enthalpy_mole(gas);
    gam_arr(j) = cp_mole(gas)/cv_mole(gas);

end

%Individual properties
% figure
% hold on
% plot(Trange,h_arr)
% plot(Trange,enthalpy(:,1),'--')
% hold off
% box off
% xlabel('Temperature [K]')
% ylabel('h [J/kg]')


%Mixture properties - Mass Units
%{
figure
hold on
plot(Trange,h_arr)
plot(Trange,sum(Y.*enthalpy,2),'--')
hold off
box off
xlabel('Temperature [K]')
ylabel('h [J/kg]')

figure
hold on
plot(Trange,gam_arr)
plot(Trange,gamma,'--')
hold off
box off
xlabel('Temperature [K]')
ylabel('$\gamma$')

figure
hold on
plot(Trange,cp_arr)
plot(Trange,cp,'--')
hold off
box off
xlabel('Temperature [K]')
ylabel('$c_p$ [J/kg-K]')

%}


%Mixture properties - Mole Units
%
figure
hold on
plot(Trange,h_arr)
plot(Trange,sum(X.*enthalpy,2),'--')
hold off
box off
xlabel('Temperature [K]')
ylabel('h [J/kmol]')

figure
hold on
plot(Trange,gam_arr)
plot(Trange,gamma,'--')
hold off
box off
xlabel('Temperature [K]')
ylabel('$\gamma$')

figure
hold on
plot(Trange,cp_arr)
plot(Trange,cp,'--')
hold off
box off
xlabel('Temperature [K]')
ylabel('$c_p$ [J/kmol-K]')

%}