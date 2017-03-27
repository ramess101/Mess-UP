clear

TC_range = [566	566.2758621	566.5517241	566.8275862	567.1034483	567.3793103	567.6551724	567.9310345	568.2068966	568.4827586	568.7586207	569.0344828	569.3103448	569.5862069	569.862069	570.137931	570.4137931	570.6896552	570.9655172	571.2413793	571.5172414	571.7931034	572.0689655	572.3448276	572.6206897	572.8965517	573.1724138	573.4482759	573.7241379	574]';

TC_plot = linspace(min(TC_range),max(TC_range),500);

PDF_TC_octane_all = dlmread('MC_TC_PDF_octane.txt');

n = length(PDF_TC_octane_all);

model = @(A,mu,sigma,x) A.*exp(-1*(x-mu).^2./(2*sigma.^2));

% SSE = @(data,mu,sigma) sum((data - normpdf(TC_range,mu,sigma)).^2);

SSE = @(data,A,mu,sigma) sum((data - model(A,mu,sigma,TC_range)).^2);

TC_sigma_guess = 0.896904208;

% figure
% hold

for i = 1:n
    
    PDF_TC = PDF_TC_octane_all(:,i);
    
    [peak, index] = max(PDF_TC);
    
    guess = [peak,TC_range(index),TC_sigma_guess];
    
    f = @(parameters) SSE(PDF_TC,parameters(1),parameters(2),parameters(3));
    
    fit = fminsearch(f,guess);
    
    TC_hat(i) = fit(2);
    TC_sigma(i) = fit(3);
    
    PDF_plot(:,i) = model(fit(1),fit(2),fit(3),TC_plot);
    
%     scatter(TC_range,PDF_TC)
%     plot(TC_plot,PDF_plot(:,i))
    
end

% hold

MC_TC_sigma_guess = mean(TC_sigma);

MC_combined = sum(PDF_plot');
MC_combined = MC_combined/sum(MC_combined);

[peak, index] = max(MC_combined);

% Fixing MC_TC to the peak

MC_TC_hat = TC_plot(index);

SSE_MC = @(A,sigma) sum((MC_combined - model(A,MC_TC_hat,sigma,TC_plot)).^2);

f = @(parameters) SSE_MC(parameters(1),parameters(2));

guess = [peak,MC_TC_sigma_guess];
   
fit = fminsearch(f,guess);

MC_TC_sigma = fit(2);

% Allowing the peak to be a parameter as well

% SSE_MC = @(A,mu,sigma) sum((MC_combined - model(A,mu,sigma,TC_plot)).^2);
% 
% f = @(parameters) SSE_MC(parameters(1),parameters(2),parameters(3));
% 
% guess = [peak,TC_plot(index),MC_TC_sigma_guess];
%    
% fit = fminsearch(f,guess);
% 
% MC_TC_hat = fit(2);
% MC_TC_sigma = fit(3);

MC_plot = model(fit(1),MC_TC_hat,MC_TC_sigma,TC_plot);

figure
hist(TC_hat)

figure
hist(TC_sigma)

figure
plot(sort(TC_sigma))

figure
hold
scatter(TC_plot,MC_combined)
plot(TC_plot,MC_plot)
    
    
