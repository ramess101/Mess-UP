function [ TC,rhoc,A,b,erra,errs ] = towhee_error_model(T,rhog,rhol,beta)
% This funciton returns the results from the towhee non-rigorous code.
% Since I am mainly just using this to find the best estimate it doesn't
% really matter that the rigorous approach was not used for the statistics.

n = 2*length(T);

rhoa = (rhog + rhol) / 2;

delrho = (rhol - rhog) / 2;
rhos = delrho;
rhos = rhos.^(1/beta);

% Since this code uses the extended length we have used a model, however, 
% I am not limited to this, if I want to I can find the standard deviations
% specific to C2, but the model might be more reliable.

% If using the model, uncomment here and below (this is without constant)
% pea = 8.12195682246548*10^-9;
% pes = 1.9346962419217932*10^-8;
% ea = 13.880645592227626;
% es = 12.18250766453814;

% If using the model, uncomment here and below (this is with constant)
b0a = 5*10^-4;
b1a = 2.25*10^-12;
b2a = 22.475;
b0s = 3.9875*10^-4;
b1s = 2.225*10^-14;
b2s = 26.55;

TC = 1.1*max(T);
TC_it = 0;

while (TC-TC_it)^2>0.0000001
    
TC_it = TC;

% If using the model, uncomment these (without constant):
% erra = pea*exp(ea*T/TC);
% errs = pes*exp(es*T/TC);

% If using the model, uncomment these (with constant):
erra = b0a + b1a*exp(b2a*T/TC);
errs = b0s + b1s*exp(b2s*T/TC);


errs = rhos * (1/beta) .* errs ./ delrho; % This should still be statistically accurate

[as, bs, sigas, sigbs] = towhee_fit(T,rhos,errs);

TC = -as/bs;
dTC = TC * sqrt((sigas/as)^2 + (sigbs/bs)^2);
a0 = as^beta;
da = ((as^beta)*beta*sigas/as);

[ar, br, sigar, sigbr] = towhee_fit(T,rhoa,erra);

b1 = -br*TC;
db = b1*sqrt((sigbr/br)^2+(dTC/TC)^2);
rhoc = ar-b1;
drhoc = sqrt(sigar^2+db^2);


end

CI = 1.3*tinv(0.95,n-4)/sqrt(n); % The 1.3 is to approximately make it the two-tailed solution in the range of 4-6 temperatures

low_TC = TC-dTC*CI;
high_TC = TC+dTC*CI;
low_rhoc = rhoc-drhoc*CI;
high_rhoc = rhoc+drhoc*CI;

A = b1/TC;
b = a0/(TC^beta);
% Chose the error model you want
% erra = pea*exp(ea*T/TC);
% errs = pes*exp(es*T/TC);

erra = b0a + b1a*exp(b2a*T/TC);
errs = b0s + b1s*exp(b2s*T/TC);

y = (rhog + rhol)/2;
z = (rhol - rhog)/2;

SSE = @(b) sum(((y - (b(2) + b(3)*(b(1)-T)))./erra).^2 + ((z - (b(4)*(b(1)-T).^beta))./errs).^2);
b_guess = [TC,rhoc,A,b];

b_fit = fminsearch(SSE,b_guess);

TC = b_fit(1);
rhoc = b_fit(2);
A = b_fit(3);
b = b_fit(4);

% figure
% hold
% plot([low_rhoc,low_rhoc],[low_TC,high_TC])
% plot([high_rhoc,high_rhoc],[low_TC,high_TC])
% plot([low_rhoc,high_rhoc],[low_TC,low_TC])
% plot([low_rhoc,high_rhoc],[high_TC,high_TC])
% scatter(rhoc,TC,'b')
% scatter(rhol,T,'g')
% scatter(rhog,T,'g')
% hold

end

