% This requires having already run "Methane_rhoL_parameter_space" (or B2)
% so that you have an eps_norm set of epsilons

STD_TC_star = 0.00081633;

TC_star_r = normrnd(TC_star,STD_TC_star,[1 length(eps_norm)]);
    
TC_PoE = eps_norm .* TC_star_r';

[TC_counts,TC_centers] = hist(TC_PoE,1000);

[TC_low,TC_high,act_TC] = integrate_histogram(TC_counts,TC_centers,0.95);

sigma_TC = ((TC_high-TC)+(TC-TC_low))/2/1.96;

PDF_TC = TC_counts/sum(TC_counts);

plot(TC_centers,PDF_TC)