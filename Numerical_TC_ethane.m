clear

TC_range = [303.75	303.8017241	303.8534483	303.9051724	303.9568966	304.0086207	304.0603448	304.112069	304.1637931	304.2155172	304.2672414	304.3189655	304.3706897	304.4224138	304.4741379	304.5258621	304.5775862	304.6293103	304.6810345	304.7327586	304.7844828	304.8362069	304.887931	304.9396552	304.9913793	305.0431034	305.0948276	305.1465517	305.1982759	305.25]';

% TC_plot = linspace(min(TC_range),max(TC_range),100);

TC_plot = linspace(303.5,305.5,500);

PDF_TC_ethane_all = dlmread('MC_TC_PDF_ethane.txt');

n = length(PDF_TC_ethane_all);

model = @(A,mu,sigma,x) A.*exp(-1*(x-mu).^2./(2*sigma.^2));

SSE = @(data,A,mu,sigma) sum((data - model(A,mu,sigma,TC_range)).^2);

TC_sigma_guess = 0.14953;

% figure
% hold

for i = 1:n
    
    PDF_TC = PDF_TC_ethane_all(:,i);
    
    [peak, index] = max(PDF_TC);
    
    guess = [peak,TC_range(index),TC_sigma_guess];
    
    f = @(parameters) SSE(PDF_TC,parameters(1),parameters(2),parameters(3));
    
    fit = fminsearch(f,guess);
    
    TC_hat(i) = fit(2);
    TC_sigma(i) = fit(3);
    
    PDF_plot(:,i) = model(fit(1),fit(2),fit(3),TC_plot);
    
%     scatter(TC_range,PDF_TC)
%     plot(TC_plot,PDF_plot(:,i))
%     
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
    
    
