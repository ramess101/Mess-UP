function [ low,high,actual_conf ] = integrate_histogram( counts,centers,confidence )
%This code attempts to integrate a histogram

dx = centers(2) - centers(1); % This code really only works when the data are evenly binned, which is normally what I do

integrand = counts*dx;

integral = sum(integrand);

norm_pdf = integrand/integral;

check_norm = sum(norm_pdf); % This should be 1

for i = 1:length(norm_pdf)
   
    norm_cdf(i) = sum(norm_pdf(1:i));
    
end

up_bound = 0.5 + confidence/2;
low_bound = 0.5 - confidence/2;

dif_up_bound = abs(norm_cdf-up_bound);
dif_low_bound = abs(norm_cdf-low_bound);

[~, I] = min(dif_up_bound);
[~, J] = min(dif_low_bound);

high = centers(I);
low = centers(J);

actual_conf = sum(norm_pdf(J:I)); % The actual confidence will not be the exact same as the desired confidence


end

