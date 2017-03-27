function [ TC,rhoc,A,b,TC_low,TC_high,rhoc_low,rhoc_high ] = towhee_rectilinear(T,rhog,rhol,errg,errl,beta )

% This is the exact same as towhee's rectilinear regression

n = 2*length(T);

rhoa = (rhog + rhol) / 2;
erra = sqrt(errg.^2 + errl.^2)/2;

delrho = (rhol - rhog) / 2;
rhos = delrho;
errs = erra;
rhos = rhos.^(1/beta);
errs = rhos * (1/beta) .* errs ./ delrho;

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

A = b1/TC;
b = a0/(TC^beta);

% CI = 1.3*tinv(0.95,n-4)/sqrt(n); % The 1.3 is to approximately make it the two-tailed solution in the range of 4-6 temperatures

TC_low = TC-dTC/sqrt(n);
TC_high = TC+dTC/sqrt(n);
rhoc_low = rhoc-drhoc/sqrt(n);
rhoc_high = rhoc+drhoc/sqrt(n);

end

