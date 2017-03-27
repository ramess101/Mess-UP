function [ PC ] = PC_from_Rackett_alt( TC, rhoc, rho_L, TR, MW)
%This function returns PC from the fits to the rectilinear and density
%scaling laws for a specific reduced temperature

D = 2/7;
Rg = 8.314472;

VC = MW/rhoc;

PC = exp(log(rhoc/rho_L)/((1-TR)^D) - log(VC/Rg/TC));

end
