function Ma = getMach_AR(gamma,AR,supersonic)
%This function computes the Mach number at the ratio A/A*

if supersonic == 1
    M = 1:0.01:10;
else
    M = 0.01:0.005:1;
end

Ma = zeros(size(AR));
for k = 1:length(AR)
    f = AR(k) - 1./M.*(2/(gamma+1)*(1+(gamma-1)/2*M.^2)).^((gamma+1)/(2*(gamma-1)));
    Ma(k) = interp1(f,M,0);

end

end
