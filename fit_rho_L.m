function [ rho_c, TC, A, B ] = fit_rho_L( T, rho_L, rho_c_g, TC_g, A_g, B_g, beta )
% This function fits just the liquid density data so that the vapor density
% data does not affect it.

rho_L_hat = @(T,b) b(1) + b(3) .* (b(2) - T) + b(4) .* (b(2) - T) .^ beta;

SE = @(b) (rho_L_hat(T,b) - rho_L).^2;

SSE = @(b) sum(SE(b));

b_fit = fminsearch(SSE,[rho_c_g, TC_g, A_g, B_g]);

rho_c = b_fit(1);
TC = b_fit(2);
A = b_fit(3);
B = b_fit(4);

end

