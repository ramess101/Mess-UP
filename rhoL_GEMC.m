function [ rhoL ] = rhoL_GEMC( T,eps,sig )
%This uses the fit of GEMC coexistence data for a rhoL curve and then
%converts the units using epsilon and sigma

T_star = T / eps; % Eps must be in units of K

rho_star = 0.31443712 + 0.34810308/2*(1.2996166-T_star) + 1.00042131/2*(1.2996166-T_star).^(0.33325329);

rhoL = rho_star / (sig^3) * 1.6605778811; % Sig must be in units of nm (the last term accounts for 1/(Avogadros * nm^3)


end

