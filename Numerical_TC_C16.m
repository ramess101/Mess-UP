clear

TC_range = [725,725.101449275362,725.202898550725,725.304347826087,725.405797101449,725.507246376812,725.608695652174,725.710144927536,725.811594202899,725.913043478261,726.014492753623,726.115942028986,726.217391304348,726.318840579710,726.420289855073,726.521739130435,726.623188405797,726.724637681159,726.826086956522,726.927536231884,727.028985507246,727.130434782609,727.231884057971,727.333333333333,727.434782608696,727.536231884058,727.637681159420,727.739130434783,727.840579710145,727.942028985507,728.043478260870,728.144927536232,728.246376811594,728.347826086957,728.449275362319,728.550724637681,728.652173913044,728.753623188406,728.855072463768,728.956521739130,729.057971014493,729.159420289855,729.260869565217,729.362318840580,729.463768115942,729.565217391304,729.666666666667,729.768115942029,729.869565217391,729.971014492754,730.072463768116,730.173913043478,730.275362318841,730.376811594203,730.478260869565,730.579710144928,730.681159420290,730.782608695652,730.884057971015,730.985507246377,731.086956521739,731.188405797102,731.289855072464,731.391304347826,731.492753623188,731.594202898551,731.695652173913,731.797101449275,731.898550724638,732]';

TC_plot = linspace(min(TC_range),max(TC_range),500);

PDF_TC_C16_all = dlmread('MC_TC_PDF_C16.txt');

n = length(PDF_TC_C16_all);

model = @(A,mu,sigma,x) A.*exp(-1*(x-mu).^2./(2*sigma.^2));

% SSE = @(data,mu,sigma) sum((data - normpdf(TC_range,mu,sigma)).^2);

SSE = @(data,A,mu,sigma) sum((data - model(A,mu,sigma,TC_range)).^2);

TC_sigma_guess = 0.849032;

% figure
% hold

for i = 1:n
    
    PDF_TC = PDF_TC_C16_all(:,i);
    
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
    
    
