% This code tries to determine what sigma and sigma2 should be for the RHS
% when using rho_L as the criteria for regression.


% First row is using the polynomial fit to the simulations, second row is
% using the error model but with a scaling factor and adjustable TC, third
% row is regressing completely new parameters for error model.

std = [0.000150006 0.000162018 0.000168543 0.000183486 0.00020333 0.000232412 0.000276759 0.000344097 0.000443843 0.000509276 0.00058711 0.000679001];
% std = [0.000194528	0.00019544	0.000196202	0.000198825	0.000204131	0.000214861	0.000236561	0.000280445	0.000369193	0.000443294	0.000548673	0.000698531];
% std = [0.000141212	0.000161368	0.000173094	0.000200449	0.000233981	0.000275085	0.00032547	0.000387232	0.000462939	0.000506981	0.000555742	0.000609728];

N = length(std);

sig = 0;
reps = 100000;

for rep = 1:reps

    sum = 0;
    
    for i = 1:N
        
        sum = sum + (normrnd(0,std(i)))^2;
        
    end
   
    sig = sig + sqrt(sum/N);
    
end

sigma = sig/reps;
sigma2 = sigma^2;
    
    