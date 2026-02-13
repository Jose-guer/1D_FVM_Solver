function [Y,dY,X,dX] = species_var(etac,phi,dx)
% This function determines the variation in species mole fraction
% and mass fraction for a given combustion efficiency curve etac(x) and
% equivalence ratio. The model is for hydrogen-air combustion and is 
% described in:
% 
% Guerrero, Jose I., and Mirko Gamba. "Quantifying combustion efficiency 
% in rotating detonation engines using MHz-rate scanned-wavelength-
% modulation spectroscopy." Combustion and Flame 285 (2026)
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

MW_H2 = 2.016;
MW_H2O = 18.0153;
MW_N2 = 28.0134;
MW_O2 = 31.999;

%Determine species mole fraction axial variation
if phi <= 1
    alpha = etac*phi;
else
    alpha = etac;
end

%From chemical balance model, determine mole fractions
Xtot = phi + (1-alpha)/2 + 1.88;
XH2O = alpha./Xtot;
XH2 = (phi-alpha)./Xtot;
XO2 = (1-alpha)/2 ./Xtot;
XN2 = 1.88./ Xtot;

%Convert to mass fractions
MWmix = XH2O.*MW_H2O + XH2.*MW_H2 + XN2.*MW_N2 + XO2.*MW_O2;
YH2O = XH2O.*MW_H2O./MWmix;
YH2 = XH2.*MW_H2./MWmix;
YN2 = XN2.*MW_N2./MWmix;
YO2 = XO2.*MW_O2./MWmix;

%Compute dXi/dx
d_XH2O = 0*etac;
d_XH2 = 0*etac;
d_XO2 = 0*etac;
d_XN2 = 0*etac;

%Compute dYi/dx
d_YH2O = 0*etac;
d_YH2 = 0*etac;
d_YO2 = 0*etac;
d_YN2 = 0*etac;

for j = 1:length(XH2O) %Mole Fractions

    if j == 1
        d_XH2O(j) = (XH2O(j+1) - XH2O(j)) /dx;
        d_XH2(j)  = (XH2(j+1)  - XH2(j)) / dx;
        d_XO2(j)  = (XO2(j+1)  - XO2(j)) / dx;
        d_XN2(j)  = (XN2(j+1)  - XN2(j)) / dx;

        d_YH2O(j) = (YH2O(j+1) - YH2O(j)) /dx;
        d_YH2(j)  = (YH2(j+1)  - YH2(j)) / dx;
        d_YO2(j)  = (YO2(j+1)  - YO2(j)) / dx;
        d_YN2(j)  = (YN2(j+1)  - YN2(j)) / dx;

    elseif j == length(XH2O)

        d_XH2O(j) = (XH2O(j) - XH2O(j-1)) /dx;
        d_XH2(j)  = (XH2(j)  - XH2(j-1)) / dx;
        d_XO2(j)  = (XO2(j)  - XO2(j-1)) / dx;
        d_XN2(j)  = (XN2(j)  - XN2(j-1)) / dx;

        d_YH2O(j) = (YH2O(j) - YH2O(j-1)) /dx;
        d_YH2(j)  = (YH2(j)  - YH2(j-1)) / dx;
        d_YO2(j)  = (YO2(j)  - YO2(j-1)) / dx;
        d_YN2(j)  = (YN2(j)  - YN2(j-1)) / dx;

    else
        d_XH2O(j) = (XH2O(j+1) - XH2O(j-1)) / (2*dx);
        d_XH2(j)  = (XH2(j+1)  - XH2(j-1)) / (2*dx);
        d_XO2(j)  = (XO2(j+1)  - XO2(j-1)) / (2*dx);
        d_XN2(j)  = (XN2(j+1)  - XN2(j-1)) / (2*dx);

        d_YH2O(j) = (YH2O(j+1) - YH2O(j-1)) / (2*dx);
        d_YH2(j)  = (YH2(j+1)  - YH2(j-1)) / (2*dx);
        d_YO2(j)  = (YO2(j+1)  - YO2(j-1)) / (2*dx);
        d_YN2(j)  = (YN2(j+1)  - YN2(j-1)) / (2*dx);

    end

end

Y = [YH2O,YH2,YN2,YO2];
dY = [d_YH2O,d_YH2,d_YN2,d_YO2];
X = [XH2O,XH2,XN2,XO2];
dX = [d_XH2O,d_XH2,d_XN2,d_XO2];

end