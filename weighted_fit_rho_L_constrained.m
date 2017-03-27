function [ rho_c, TC, A, B ] = weighted_fit_rho_L_constrained( T, rho_L, rho_c_g, TC_g, A_g, B_g, beta,TC_low )
% This function fits just the liquid density data so that the vapor density
% data does not affect it.

sigma_L = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(T/TC_g));

rho_L_hat = @(T,b) b(1) + b(3) .* (b(2) - T) + b(4) .* (b(2) - T) .^ beta;

WSE = @(b) ((rho_L_hat(T,b) - rho_L)./sigma_L).^2;

WSSE = @(b) sum(WSE(b));

if TC_g < TC_low
    
    TC_g = TC_low;
    
end

b_g = [rho_c_g, TC_g, A_g, B_g];

b_fit = fmincon(WSSE,b_g,[],[],[],[],[0,TC_low,0,0],[]);

rho_c = b_fit(1);
TC = b_fit(2);
A = b_fit(3);
B = b_fit(4);

end

