function [ rhoc, TC, A, B, beta ] = fit_rho_L_beta( T, rho_L, rhoc_g, TC_g, A_g, B_g, beta_g )
% This function fits just the liquid density data so that the vapor density
% data does not affect it.

rho_L_hat = @(T,b) b(1) + b(3) .* (b(2) - T) + b(4) .* (b(2) - T) .^ b(5);

SE = @(b) (rho_L_hat(T,b) - rho_L).^2;

SSE = @(b) sum(SE(b));

b_fit = fminsearch(SSE,[rhoc_g, TC_g, A_g, B_g, beta_g]);

rhoc = b_fit(1);
TC = b_fit(2);
A = b_fit(3);
B = b_fit(4);
beta = b_fit(5);

end

