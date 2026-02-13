% This script implements a first-order accurate, Finite Volume Method (FVM) 
% solver with a Roe flux for the quasi-1D compressible flow equations.
% Becuase the Roe flux is an approximate solution to Riemann problem, the
% solver is shock capturing.
% 
% This particular implementation allows for area changes, the effects of 
% friction, and chemical reactions. As an example, a combustor with a 
% converging-divering nozzle is hard coded.
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

clc
clear all
close all

% Solution settings
CFL  = 0.5;                % courant number
Max_iterations = 20000;
resmax = 10^-4;            % maximum error

cf = 0;                  % skin friction coefficient

% Combustion efficiency varies linearly from 0 at x0 to etac_f at xf
% x0 < 0, combustor section from x=x0 to x=xf=0
etac_f = 1;
x0 = -0.5;
xf = 0;

%---------------------------------------
%H2-Air Adiabatic Flame Temp.

% phi = 0.6;
% Teq = 2000; %phi = 0.6

% phi = 0.8;
% Teq = 2341; %phi = 0.8

phi = 1;
Teq = 2608; %phi = 1.0

% phi = 1.2;
% Teq = 2559.02; %phi = 1.2

%---------------------------------------
% Boundary Conditions

% Subsonic inflow
Ts = 500;      %static temperature
Ps = 80 * 101325;

% Supersonic outflow
exit_BC_supersonic = 1;

% Subsonic outflow
% exit_BC_supersonic = 0;
% Pe = 101325;

%---------------------------------------
% Nozzle profile

Lc = x0;             % Combustor Inlet
L1       = 1;        % Length converging section
L2       = 2;        % Length diverging section
Ae_to_Ath = 8.1684;  % throat area / exit area
D_inlet  = 0.75;     % m
D_exit   = 0.5;      % m

x = linspace(Lc, L1 + L2, 500);
x = x';
dx = x(2) - x(1);

ind = x >= 0; % Nozzle
xc = x(~ind); %Combustor
D = laval_nozzle_C1(x(ind), L1, L2, Ae_to_Ath, D_inlet, D_exit);
D0 = D(1);
D = [D0*ones(size(xc)); D];

A = pi/4*D.^2; %Area, m^2
A_inlet = A(1);
A_exit = A(end);
A_th = min(A);

%---------------------------------------

LnA = log(A);
dLnAdx = 0*x;
k = length(x);
for i = 1:k
    if i == 1
        dLnAdx(i) = (LnA(i+1) - LnA(i)) / dx;
    elseif i == k
        dLnAdx(i) = (LnA(i) - LnA(i-1)) / dx;
    else
        dLnAdx(i) = (LnA(i+1) - LnA(i-1)) / (2*dx);
    end
end

F = 4*cf./D;

figure
hold on
plot(x,D/2,'k')
plot(x,-D/2,'k')
axis equal
box off
xlabel('x')
ylabel('Radius')

%% Initalize Flowfield

ind = x < 0; %combustor region

%--------------------------------------------------------------------------
% Pressure 

P = 0 * x;
P(ind) = Ps; %constant pressure combustion

%Estimate exit pressure
Ma = getMach_AR(1.3,Ae_to_Ath,1);
g = 1.4;
Pe = Ps / (1 + (g-1)/2 * Ma^2)^(g/(g-1));

N = sum(~ind);
P(~ind) = linspace(Ps,Pe,N); %Pressure varries linearly from combustor to exit

%--------------------------------------------------------------------------
% Temperature 

T = 0 * x;

%Temperature varries linearly from the Ts at the inlet to the T_adiab. at the
%combustor outlet
N = sum(ind);
T(ind) = linspace(Ts,Teq,N);

% Option #1 - Assume temperature is constant afterwards
T(~ind) = Teq;

% Option #2 - Assume isentropic flow
% T(~ind) = Teq .* (P(~ind)/Ps).^ ((1.4-1)/1.4);
%}

%--------------------------------------------------------------------------
% Flowrate  

M_in = 0.01;
V_in = M_in*sqrt(1.4*287*Ts);
mdot = Ps/(287*Ts) * V_in * A_inlet;

%**************************************************************************
%                       Etac and Thermoprop
%**************************************************************************

etac = etac_f *(x-x0)/(xf-x0);
ind = x > xf;
etac(ind) = etac_f;

[Y,dY,~,~] = species_var(etac,phi,dx);

mol_units = 0;

if etac_f ~= 0
    [gamma, enthalpy, R] = thermdat_H2Air(T,Y,mol_units);
else
    gamma = 1.4 * ones(size(x));
    R = 287 * ones(size(x));
end

rho = P ./ (R.*T);
V = mdot./(rho.*A);
E = R.*T./(gamma-1) + 0.5.*V.^2;

%**************************************************************************
%           Initial condition of solution vectors
%**************************************************************************

U1 = rho.*A;
U2 = rho.*V.*A;
U3 = rho.*E.*A;

U = [U1, U2, U3];

figure; set(gcf,'Position',[293 149 1223 791])
subplot(3,1,1)
plot(x,U1,'r')
box off
ylabel('$U_1 = \rho A$')
title('Initial Condition')
subplot(3,1,2)
plot(x,U2,'b')
box off
ylabel('$U_2 = \dot{m}$')
subplot(3,1,3)
plot(x,U3,'k')
box off
xlabel('x (m)')
ylabel('$U_3 = \rho E A$')

figure; set(gcf,'Position',[293 149 1223 791])
subplot(3,1,1)
plot(x,rho,'r')
box off
ylabel('$\rho$')
title('Initial Condition')
subplot(3,1,2)
plot(x,T,'b')
box off
ylabel('$T$')
subplot(3,1,3)
plot(x,P,'k')
box off
xlabel('x (m)')
ylabel('$P$')

% return

%% 1D-FVM

%
res = 1;
t = 0;
nstep = 0;

F_minus_half = 0*U;
F_plus_half = 0*U;

J = 0*U;
%}

tic

iter = 1;

while res > resmax

    %**********************************************************************
    %                        Time step
    %**********************************************************************

    P = (gamma-1)./A .*( U(:,3) - 0.5 .* U(:,2).^2 ./ U(:,1) );
    ind = P < 0; P(ind) = eps;
    rho = U(:,1)./A;
    a = sqrt(gamma.*P./rho);
    u = U(:,2)./U(:,1);

    dta = (CFL*dx)./(abs(u)+a);
    dt = min(dta);
    
    t = t + dt;
    nstep = nstep+1;


    %**********************************************************************
    %                        Roe Flux
    %**********************************************************************

    for i = 2:k-1 %interior points
        F_plus_half(i,:)  = ROEflux2(U(i,:), U(i+1,:), A(i), A(i+1), gamma(i), gamma(i+1));
        F_minus_half(i,:) = ROEflux2(U(i-1,:), U(i,:),A(i-1), A(i), gamma(i-1), gamma(i));
    end

    %**********************************************************************
    %                         Source Terms
    %**********************************************************************

    J(:,2) = (gamma-1).* ( U(:,3) - 0.5 .*  U(:,2).^2 ./ U(:,1) ) .* dLnAdx - 0.5.* U(:,2).^2 ./ U(:,1) .* F;

    % *********** Heat addition ***********************

    if etac_f ~= 0

        qdot = - 1 .* enthalpy .* U(:,2).* dY;
        qdot = sum(qdot,2);

        %include frictional terms in energy equation
        % J(:,3) =  qdot - 0.5.* U(:,2).^3 ./ U(:,1).^2 .* F;

        %exclude fricitonal terms in energy equation
        J(:,3) =  qdot;

    else
        %include frictional terms in energy equation
        % J(:,3) = - 0.5.* U(:,2).^3 ./ U(:,1).^2 .* F;

        %exclude fricitonal terms in energy equation
        J(:,3) =  0;
    end
    % ************** End heat addition  *****************


    U(:,1) = U(:,1) - dt/dx .* (F_plus_half(:,1) - F_minus_half(:,1)) + dt .* J(:,1);
    U(:,2) = U(:,2) - dt/dx .* (F_plus_half(:,2) - F_minus_half(:,2)) + dt .* J(:,2);
    U(:,3) = U(:,3) - dt/dx .* (F_plus_half(:,3) - F_minus_half(:,3)) + dt .* J(:,3);

    R2 = - dt/dx .* (F_plus_half(:,2) - F_minus_half(:,2)) + dt .* J(:,2);
    eps = 10^-10;

    % res = max(abs( R2./(U(:,2)+eps) ));
    res = max(abs( R2 ));

    res_mom(iter) = res;
    iter = iter + 1;

    %**********************************************************************
    %                        Boundary Conditions
    %**********************************************************************

    % *************** Subsonic Inlet ***************

    % *** Option 1, Fix Ps, Ts ***
    rhos = Ps./(R(1).*Ts);
    U(1,1) = rhos(1)*A(1);
    U(1,2) = 2*U(2,2)-U(3,2); %U2
    V_inlet = U(1,2)/U(1,1);
    E_inlet = R(1).*Ts./(gamma(1)-1) + 0.5 * V_inlet.^2;
    U(1,3) = U(1,1).*E_inlet;

    % *** Option 2, Fix mdot, Ts ***
    % U(1,1) = 2*U(2,1)-U(3,1); %U1
    % U(1,2) = mdot;
    % V_inlet = U(1,2)/U(1,1);
    % E_inlet = R*Ts/(gamma-1) + 0.5 * V_inlet.^2;
    % U(1,3) = U(1,1).*E_inlet;

    % *** Option 3, Fix mdot, Ps ***
    % U(1,1) = 2*U(2,1)-U(3,1); %U1
    % U(1,2) = mdot;
    % U(1,3) = Ps*A(1)/(gamma-1) + 0.5*U(1,2)^2/U(1,1);

    if exit_BC_supersonic == 0
        % *************** Subsonic Outlet ***************
        U(k,1) = 2*U(k-1,1)-U(k-2,1);
        U(k,2) = 2*U(k-1,2)-U(k-2,2);
        U(k,3) = Pe*A(k)./(gamma(k)-1) + 0.5 * U(k,2).^2 ./ U(k,1);

    else
        % *************** Supersonic Outlet ***************
        U(k,1) = 2*U(k-1,1)-U(k-2,1);
        U(k,2) = 2*U(k-1,2)-U(k-2,2);
        U(k,3) = 2*U(k-1,3)-U(k-2,3);
    end

    %**********************************************************************
    %                        Species variation
    %**********************************************************************

    %use old gamma values to compute a new temperature
    if etac_f ~= 0
        T = (gamma-1)./R .* (U(:,3)./U(:,1) - 0.5.*(U(:,2)./U(:,1)).^2 );
        [gamma, enthalpy, R] = thermdat_H2Air(T,Y,mol_units);
        [~,ind] = min(A);
        gam_th = gamma(ind);
        gamma(ind:end) = gam_th;
    end

    if nstep == Max_iterations
        break;
    end

end

toc

% return

%% Post Processing

U1 = U(:,1); U2 = U(:,2); U3 = U(:,3);


rho = U1./A;
V = U2./U1;
E = U3./U1;
T = (E - 0.5*V.^2) .* (gamma - 1) ./ R;
P = rho .* R .*T;
a = sqrt(gamma.*R.*T);
M = V./a;
mf = U2;

Tt = T.*(1 + (gamma-1)/2 .* M.^2);
Pt = P.*(Tt./T).^(gamma./(gamma-1));

[~, enthalpy, ~] = thermdat_H2Air(T,Y,mol_units);

qdot = - 1 .* enthalpy .* mf.* dY ;
qdot = sum(qdot,2);

h0 = sum(Y.*enthalpy,2) + 0.5 * V.^2;


%%

figure
semilogy(res_mom)
xlabel('iteration')
ylabel('log(Residual)')
box off

%% Plots

m = 100;

figure
hold on
plot(x,P/101325)
hold off
xlabel('x [m]')
ylabel('P [Pa]')

figure
hold on
plot(x,T)
hold off
xlabel('x [m]')
ylabel('T [K]')

figure
hold on
plot(x,rho)
hold off
xlabel('x [m]')
ylabel('$\rho$ [kg/m$^3$]')

figure
hold on
plot(x,M)
hold off
xlabel('x [m]')
ylabel('Mach number')

figure
hold on
plot(x,mf)
xlabel('x [m]')
ylabel('$\dot{m}$ [kg/s]')
box off
ylim([50,150])
xlim([-0.5,3])

figure
hold on
plot(x,gamma)
xlabel('x [m]')
ylabel('$\gamma$')
box off

%% Surf Plots

n = 25;

% Grid distribution in the nozzle
for i = 1:k
    y = linspace(-D(i)/2,D(i)/2,n);
    for j = 1:n
        xx(i,j) = x(i);
        yy(i,j) = y(j);
    end
end

for i = 1:k
    for j = 1:n
        rho1(i,j) = rho(i);
        T1(i,j) = T(i);
        P1(i,j) = P(i);
        M1(i,j) = M(i);
        V1(i,j) = V(i);
    end
end

figure
surf(xx,yy,M1);view(2)
shading interp
colormap(jet)
colorbar
% axis off
grid off
title('Mach Number')
xlabel('x [m]')
ylabel('r [m]')
axis equal

figure
surf(xx,yy,T1);view(2)
shading interp
colormap(jet)
c=colorbar;
% axis off
grid off
title('Temperature')
ylabel('r [m]')
xlabel('x [m]')
axis equal

figure
surf(xx,yy,P1/101325);view(2)
shading interp
colormap(jet)
colorbar
% axis off
grid off
title('Pressure')
xlabel('x [m]')
ylabel('r [m]')
axis equal
