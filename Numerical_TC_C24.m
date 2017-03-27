clear

TC_range = [800,800.769230769231,801.538461538462,802.307692307692,803.076923076923,803.846153846154,804.615384615385,805.384615384615,806.153846153846,806.923076923077,807.692307692308,808.461538461539,809.230769230769,810,810.769230769231,811.538461538462,812.307692307692,813.076923076923,813.846153846154,814.615384615385,815.384615384615,816.153846153846,816.923076923077,817.692307692308,818.461538461539,819.230769230769,820,820.769230769231,821.538461538462,822.307692307692,823.076923076923,823.846153846154,824.615384615385,825.384615384615,826.153846153846,826.923076923077,827.692307692308,828.461538461539,829.230769230769,830]';

TC_plot = linspace(min(TC_range),max(TC_range),500);

PDF_TC_C24_all = dlmread('MC_TC_PDF_C24.txt');

n = length(PDF_TC_C24_all);

model = @(A,mu,sigma,x) A.*exp(-1*(x-mu).^2./(2*sigma.^2));

% SSE = @(data,mu,sigma) sum((data - normpdf(TC_range,mu,sigma)).^2);

SSE = @(data,A,mu,sigma) sum((data - model(A,mu,sigma,TC_range)).^2);

TC_sigma_guess = 2;

figure
hold

for i = 1:n
    
    if i == 24
        
        pause(0)
        
    end
    
    PDF_TC = PDF_TC_C24_all(:,i);
    
    [peak, index] = max(PDF_TC);
    
    guess = [peak,TC_range(index),TC_sigma_guess];
    
    f = @(parameters) SSE(PDF_TC,parameters(1),parameters(2),parameters(3));
    
    fit = fminsearch(f,guess);
    
    TC_hat(i) = fit(2);
    TC_sigma(i) = fit(3);
    
    PDF_plot(:,i) = model(fit(1),fit(2),fit(3),TC_plot);
    
    scatter(TC_range,PDF_TC)
    plot(TC_plot,PDF_plot(:,i))
     
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
    
mean(TC_sigma)
    
