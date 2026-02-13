function D = laval_nozzle_C1(x, L1, L2, Ae_to_Ath, D_inlet, D_exit)
%LAVAL_NOZZLE_C1  Bell-shaped converging-diverging nozle profile.
%   D = LAVAL_NOZZLE_C1(x, L1, L2, Ae_to_Ath, D_inlet, D_exit)
%   returns the diameter D(x)
%
%   Axial coordinate convention:
%       x = 0        : nozzle inlet
%       x = L1       : throat
%       x = L1 + L2  : nozzle exit
%
%   INPUTS:
%       x        : vector of axial coordinates [m]
%       L1       : inlet-to-throat length [m]
%       L2       : throat-to-exit length [m]
%       Ae_to_Ath : exit-to-throat area ratio
%       D_inlet  : inlet diameter [m]
%       D_exit   : exit diameter [m]
%
%   OUTPUT:
%       D        : diameter at each x [m]
%
%   The contour is C1 continuous with:
%   dD/dx = 0 at inlet, throat, and exit

    % --- throat diameter from area ratio ---
    D_throat = D_exit * sqrt(1/Ae_to_Ath);

    % --- allocate output ---
    D = zeros(size(x));

    % --- converging section: 0 <= x <= L1
    % Hermite blend: D(s) = Dinlet + (Dthroat - Dinlet)*(3s^2 - 2s^3)
    % with s = x/L1, giving dD/dx = 0 at s=0 (inlet) and s=1 (throat).
    idx1 = (x <= L1);
    if any(idx1)
        s = x(idx1) ./ L1;           % 0 -> 1
        H = 3.*s.^2 - 2.*s.^3;       % smooth step
        D(idx1) = D_inlet + (D_throat - D_inlet) .* H;
    end

    % --- diverging (bell) section: L1 < x <= L1+L2
    % Same Hermite blend from throat to exit:
    % D(t) = Dthroat + (Dexit - Dthroat)*(3t^2 - 2t^3)
    % with t = (x - L1)/L2, giving dD/dx = 0 at throat and exit.
    idx2 = (x > L1);
    if any(idx2)
        t = (x(idx2) - L1) ./ L2;    % 0 -> 1
        H = 3.*t.^2 - 2.*t.^3;
        D(idx2) = D_throat + (D_exit - D_throat) .* H;
    end
end
