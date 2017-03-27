function [ rhoc, ZC, TC, D ] = weighted_fit_rho_L_Rackett( T, rho_L, rhoc_g, TC_g, PC_g, MW )
% This function fits just the liquid density data so that the vapor density
% data does not affect it.

ZC_g = PC_g * MW / rhoc_g / 8.314472 / TC_g;

sigma_L = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(T/TC_g));

rho_L_hat = @(T,b) b(1) .* b(2) .^ (-(1-T./b(3)).^(b(4)));

WSE = @(b) ((rho_L_hat(T,b) - rho_L)./sigma_L).^2;

WSSE = @(b) sum(WSE(b));

b_fit = fminsearch(WSSE,[rhoc_g, ZC_g, TC_g, 2/7]);

rhoc = b_fit(1);
ZC = b_fit(2);
TC = b_fit(3);
D = b_fit(4);

end

