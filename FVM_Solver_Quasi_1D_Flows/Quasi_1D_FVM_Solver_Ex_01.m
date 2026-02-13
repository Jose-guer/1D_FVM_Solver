% This script implements a first-order accurate, Finite Volume Method (FVM) 
% solver with a Roe flux for the quasi-1D compressible flow equations.
% Becuase the Roe flux is an approximate solution to Riemann problem, the
% solver is shock capturing.
% 
% This particular implementation allows for area changes and the effects of 
% friction. Four examples are hard coded and the FVM solution is compared 
% to the exact solution. 
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

clc
clear all
close all

%Solution settings
dx = 0.01;              % grid size
CFL  = 0.5;              % courant number
Max_iterations = 10000;
resmax = 10^-4;          % maximum error

%Gas properties
MW = 28.95;            % Molecular weight
gamma = 1.4;           % Ratio of specific heats
R =  8314.5/MW;        % gas constant, J/kg-K

%----- Boundary Conditions -------

% Example 01 - Isentropic nozzle
%{
cf = 0;                % skin friction coefficient
load('Isentropic_Nozzle.mat')

% Subsonic inflow
Ts = 1986.87;      %static temperature
Ps = 7921305.147;
mdot = 535.1376;

% Subsonic outflow
Pe = 41348.35;
exit_BC_supersonic = 1;
%}

% Example 02 - Nozzle with friction
%{
cf = 0.1;                % skin friction coefficient
load('Nozzle_w_Friction.mat')

% Subsonic inflow
Ts = 1990.9831;      %static temperature
Ps = 7978811.211;
mdot = 445.8079;

% Subsonic outflow
Pe = 133855.839;
exit_BC_supersonic = 1;
%}

% Example 03 -  Shock x = 2m, no friction
%{
cf = 0;                % skin friction coefficient
load('Nozzle_w_Shock.mat')

% Subsonic inflow
Ts = 551.357;      %static temperature
Ps = 673766.319;
mdot = 86.4065;

% Subsonic outflow
Pe = 228880.547;
exit_BC_supersonic = 0;
%}

% Example 04 -  Shock x = 2m, with cf = 0.1
%
cf = 0.1;                % skin friction coefficient
load('Nozzle_w_Shock_and_Friction.mat')

% Subsonic inflow
Ts = 552.50;      %static temperature
Ps = 678657.64;
mdot = 71.983;

% Subsonic outflow
Pe = 146201.46;
exit_BC_supersonic = 0;
%}

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discretization of distance along the nozzle
x = 0:dx:3;
x = x';
k = length(x);

% Nozzle profile
D = sqrt(0.093 + 0.25*(x-0.9144).^2);
A = pi/4*D.^2; %Area, m^2
A_inlet = A(1);
A_th = min(A);
LnA = log(A);

dLnAdx = 0*x;

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
xlabel('x')
ylabel('Diameter')

% return

%% Initalize Flowfield

%**************************************************************************
%Exact Solution
%**************************************************************************

xq = dat.x; %m
rho_q = dat.rho; %kg/m^3
T_q = dat.T; %K
P_q = dat.P;
M_q = dat.M;
mf_0 = dat.mdot; %kg/s
% Pe = dat.Pe; %Pa

%**************************************************************************
% Initalize Flowfield
%**************************************************************************

N = length(x);

%Option 1: assume a constant density everywhere
%
rhox = Ps/(R*Ts)*ones(N,1);
Vx = mdot./(rhox.*A);
Ex = R.*Ts./(gamma-1) + 0.5.*Vx.^2;
%}

%Option 2: assume an isentropic solution
%{
AR = A./A_th;
Ma = getMach_Ar(gamma,AR,0);

% ind = find(A == A_th);
% Ma_sub = getMach_Ar(gamma,AR(1:ind-1),0);
% Ma_sup = getMach_Ar(gamma,AR(ind:end),0);
% Ma = [Ma_sub;Ma_sup];

M_in = Ma(1);

T0 = Ts * (1 + (gamma-1)/2 *M_in.^2);
P0 = Ps * (T0./Ts).^(gamma/(gamma-1));

Tx = T0 ./ (1 + (gamma-1)/2 * Ma.^2);
Px = P0 .* (Tx./T0).^(gamma/(gamma-1));
rhox = Px./(R.*Tx);

ax = sqrt(gamma*R*Tx);
Vx = Ma.*ax;
Ex = R.*Tx./(gamma-1) + 0.5*Vx.^2;
%}

%Option 3: use exact solution
%{
rhox = interp1(xq,rho_q,x,"linear","extrap");
Tx = interp1(xq,T_q,x,"linear","extrap");
Vx = mdot./(rhox.*A);
Ex = R.*Tx./(gamma-1) + 0.5.*Vx.^2;
%}

U1 = rhox.*A;
U2 = rhox.*Vx.*A;
U3 = rhox.*Ex.*A;

U = [U1, U2, U3];

% figure; set(gcf,'Position',[293 149 1223 791])
figure;
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

% return

%% 1D-FVM

res = 1;
t = 0;
nstep = 0;

F_minus_half = 0*U;
F_plus_half = 0*U;

J = 0*U;

while res > resmax

    %**********************************************************************
    %                        Time step
    %**********************************************************************

    P = (gamma-1)./A .*( U(:,3) - 0.5 .* U(:,2).^2 ./ U(:,1) );
    rho = U(:,1)./A;
    a = sqrt(gamma.*P./rho);
    u = U(:,2)./U(:,1);

    dta = CFL*dx./(abs(u)+a);
    dt = min(dta);

    nstep = nstep+1;

    %**********************************************************************
    %                        Roe Flux
    %**********************************************************************

    for i = 2:k-1 %interior points
        F_plus_half(i,:)  = ROEflux(U(i,:),U(i+1,:),A(i),A(i+1),gamma) ;
        F_minus_half(i,:) = ROEflux(U(i-1,:),U(i,:),A(i-1),A(i),gamma) ;
    end

    J(:,2) = (gamma-1).* ( U(:,3) - 0.5 .*  U(:,2).^2 ./ U(:,1) ) .* dLnAdx - 0.5.* U(:,2).^2 ./ U(:,1) .* F;

    U(:,1) = U(:,1) - dt/dx .* (F_plus_half(:,1) - F_minus_half(:,1)) + dt .* J(:,1);
    U(:,2) = U(:,2) - dt/dx .* (F_plus_half(:,2) - F_minus_half(:,2)) + dt .* J(:,2);
    U(:,3) = U(:,3) - dt/dx .* (F_plus_half(:,3) - F_minus_half(:,3)) + dt .* J(:,3);

    R2 = - dt/dx .* (F_plus_half(:,2) - F_minus_half(:,2)) + dt .* J(:,2);
    res = max(abs( R2 ));

    %**********************************************************************
    %                        Boundary Conditions
    %**********************************************************************

    % ********** Subsonic Inlet **********

    % *** Option 1, Fix Ps, Ts ***
    rhos = Ps/(R*Ts);
    U(1,1) = rhos*A(1);
    U(1,2) = 2*U(2,2)-U(3,2); %U2
    V_inlet = U(1,2)/U(1,1);
    E_inlet = R*Ts/(gamma-1) + 0.5 * V_inlet.^2;
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
        % ********** Subsonic Outlet **********
        U(k,1) = 2*U(k-1,1)-U(k-2,1);
        U(k,2) = 2*U(k-1,2)-U(k-2,2);
        U(k,3) = Pe*A(k)/(gamma-1) + 0.5 * U(k,2).^2 ./ U(k,1);
    else
        % ********** Supersonic Outlet **********
        U(k,1) = 2*U(k-1,1)-U(k-2,1);
        U(k,2) = 2*U(k-1,2)-U(k-2,2);
        U(k,3) = 2*U(k-1,3)-U(k-2,3);
    end

    if nstep == Max_iterations
        break;
    end

end

%% Post Processing

U1 = U(:,1); U2 = U(:,2); U3 = U(:,3);

rho = U1./A;
V = U2./U1;
E = U3./U1;
T = (E - 0.5*V.^2) * (gamma - 1) / R;
P = rho .* R .*T;
a = sqrt(gamma.*R.*T);
M = V./a;
% mf = rho.*A.*V;
% mf = P.*A.*M .* sqrt(gamma./(R.*T));
mf = U2;

Tt = T.*(1 + (gamma-1)/2 * M.^2);
Pt = P.*(Tt./T).^(gamma/(gamma-1));

%% Plots

m = 20;

figure
hold on
plot(x,P/101325)
plot(xq(1:m:end),P_q(1:m:end)/101325,'*','MarkerSize',6)
hold off
legend('1D-FVM','Exact','box','on')
xlabel('x [m]')
ylabel('P [atm]')


figure
hold on
plot(x,T)
plot(xq(1:m:end),T_q(1:m:end),'*','MarkerSize',6)
hold off
legend('1D-FVM','Exact','box','on')
xlabel('x [m]')
ylabel('T [K]')

figure
hold on
plot(x,M)
plot(xq(1:m:end),M_q(1:m:end),'*','MarkerSize',6)
hold off
legend('1D-FVM','Exact','box','on')
xlabel('x [m]')
ylabel('Mach number')
legend('1D-FVM','Exact','box','on','location','northwest')

figure
hold on
plot(x,mf)
yline(mf_0,'k--','LineWidth',2)
xlabel('x [m]')
ylabel('$\dot{m}$ [kg/s]')
box off
legend('1D-FVM','Exact','box','on')
ylim([mf_0-200 ,mf_0 + 200])


%% Surf Plots
%

n = 25;

% Grid distribution in the nozzle
for i = 1:k
    y = linspace(-A(i)/2,A(i)/2,n);
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
colormap(jet(256))
colorbar
axis off
drawnow
title('Mach number contour')
%}
