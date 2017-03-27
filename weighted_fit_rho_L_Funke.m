function [ b_fit, WSSE_fit ] = weighted_fit_rho_L_Funke( T, rho_L, rho_c_g, TC_g )
% This function fits just the liquid density data so that the vapor density
% data does not affect it.

sigma_L = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(T/TC_g));

rho_L_hat = @(T,b) b(1) .* exp(b(3) .* (1-T./b(2)).^(0.329) + b(4) .* (1-T./b(2)).^(4/6) + b(5) .* (1-T./b(2)).^(8/6) + b(6) .* (1-T./b(2)).^(19/6)); 

WSE = @(b) ((rho_L_hat(T,b) - rho_L)./sigma_L).^2;

WSSE = @(b) sum(WSE(b));

Y = (log(rho_L./rho_c_g))';

x = (1-T./TC_g)';

X = [x.^(0.329), x.^(4/6), x.^(8/6), x.^(19/6)]; 

b_g = (X'*X)\(X'*Y);

b_fit = fminsearch(WSSE,[rho_c_g, TC_g, b_g(1), b_g(2), b_g(3), b_g(4)]);

WSSE_fit = WSSE(b_fit);

end

