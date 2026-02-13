function Roe = ROEflux2(UL,UR,AL,AR,gamL,gamR)
% This function computes the Roe flux for chemically reacting flows, or
% flows in which gamma is not assumed to be constant.
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

% Left state
rhoL = UL(1)/AL;
uL = UL(2)./UL(1);
EL = UL(3)./UL(1);
pL = (gamL-1)/AL * ( UL(3) - 0.5*UL(2)^2/UL(1) );
HL = EL + pL/rhoL;
aL = (gamL-1)*(HL- 0.5*uL^2); %sound speed, squared

% Right state
rhoR = UR(1)/AR;
uR = UR(2)./UR(1);
ER = UR(3)./UR(1);
pR = (gamR-1)/AR * ( UR(3) - 0.5*UR(2)^2/UR(1) );
HR = ER + pR/rhoR;
aR = (gamR-1)*(HR- 0.5*uR^2); %sound speed, squared

% First compute the Roe Averages
r = sqrt(rhoR/rhoL);
u = (uL+r*uR)/(1+r);
H = (HL+r*HR)/(1+r);
rho = sqrt(rhoR*rhoL);

if gamL == gamR
    gamma = gamL;
    a = sqrt( (gamma-1)*(H-u*u/2) );
else
    a = sqrt((aL+r*aR)/(1+r));
end

% Differences in primitive variables
dr = rhoR - rhoL;
du = uR - uL;
dP = pR - pL;

% Wave strength (Characteristic Variables).
alpha = [(dP-rho*a*du)/(2*a^2); dr-dP/(a^2); (dP+rho*a*du)/(2*a^2)];

% Wave speeds (Eigenvalues)
ws =  [u-a; u; u+a];

% Harten's Entropy Fix
eps = 0.05*a;
for j = 1:3
    if abs(ws(j)) < eps
        ws(j) = (ws(j)^2 + eps^2)/(2*eps);
    end
end

%Absolute value of Eigenvalues
ws = abs(ws);

% Right Eigenvectors
R = [  1  ,  1  ,  1  ;
    u-a ,  u  , u+a ;
    H-u*a,u^2/2,H+u*a];

% Compute left and right fluxes using primitive variables
FL=[rhoL.*uL; rhoL.*uL.^2+pL; uL.*(rhoL.*EL+pL)];
FR=[rhoR.*uR; rhoR.*uR.^2+pR; uR.*(rhoR.*ER+pR)];

Aint = 0.5*(AL + AR); %interface area

% Roe flux F^_(i+1/2), use matrix multiplication
Roe = Aint * ( FL + FR  -  R*(ws.*alpha) )/2;
Roe = Roe';

end