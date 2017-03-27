[TC_fit, pc_fit, A_fit, b_fit, sy, sz] = towhee_error_model(T,pv,pL,beta);

n = 2*length(T);

p = 4;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sy).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sz).^2; % This includes uncertainties

SSE_fit = sum(SE_fit);

sigma = SSE_fit/(n-p);

RHS = sigma * (n + p * (finv(0.95,p,n-p)-1)); % No 0.95^p because I want the true 95%

% A_range = linspace(0.9999,1.0001,30)*A_fit;
% b_range = linspace(0.9999,1.0001,30)*b_fit;
% pc_range = linspace(0.9999,1.0001,30)*pc_fit;
% TC_range = linspace(0.9999,1.001,30)*TC_fit;

% For the MCO runs

A_range = linspace(0.95,1.05,20)*A_fit;
b_range = linspace(0.99,1.01,20)*b_fit;
pc_range = linspace(0.98,1.02,20)*pc_fit;
TC_range = linspace(0.997,1.003,20)*TC_fit;

s=1;
for g=1:length(A_range)
for h=1:length(b_range)
for i=1:length(pc_range)
for j=1:length(TC_range)
SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
SSE = sum(SE);
if SSE < RHS
A_ext(s) = A_range(g);
b_ext(s) = b_range(h);
pc_ext(s) = pc_range(i);
TC_ext(s) = TC_range(j);
s=s+1;
if SSE < SSE_fit % This is to make sure that we actually have the global minimum. If not, rerun the analysis with the new optimum (after plugging in as a new Mathcad guess)
new_best_fit = [A_range(g) b_range(h) pc_range(i) TC_range(j)];
SSE_fit = SSE;
end
end
end
end
end
end
% This eliminates the superfluous zero elements
A_ext = A_ext(:,1:(s-1));
b_ext = b_ext(:,1:(s-1));
pc_ext = pc_ext(:,1:(s-1));
TC_ext = TC_ext(:,1:(s-1));
A_low_temp = min(A_ext);
A_high_temp = max(A_ext);
b_low_temp = min(b_ext);
b_high_temp = max(b_ext);
pc_low_temp = min(pc_ext);
pc_high_temp = max(pc_ext);
TC_low_temp = min(TC_ext);
TC_high_temp = max(TC_ext);

A_low_temp/A_fit
A_high_temp/A_fit
b_low_temp/b_fit
b_high_temp/b_fit
pc_low_temp/pc_fit
pc_high_temp/pc_fit
TC_low_temp/TC_fit
TC_high_temp/TC_fit