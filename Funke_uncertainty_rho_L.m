clear
%This script uses all the rho_L data for ethane and fits the rectilinear +
%scaling form only to rho_L to find the uncertainty in rho_L

Temp_exp = [88.71 90.37 91 94.26 95 99.82 100 105.37 110 110.93 116.48 120 122.04 127.59 130 133.15 140 150 160 170 175 185 195 200 210 220 230 240 250 260 265 270 275 280 283.2 283.2 285 288.19 288.19 290 293.18 293.18 295 295.67 295.67 298 298.17 298.17 300 300.66 300.66 301 302 302.16 302.16 303 303.16 304.15 304.65 305.15 303.5 304 304.4 304.7 304.9 305 305.1 305.2 305.25 305.275]; 
DIPPR_exp = [0.6519176 0.65152669 0.650847108 0.64611409 0.646468916 0.64040079 0.640987155 0.63450707 0.630014612 0.62882384 0.62290005 0.618981929 0.61700633 0.6110224 0.60786505 0.60509861 0.59661887 0.585207305 0.573594271 0.561722635 0.555675558 0.543319795 0.530558087 0.52399982 0.51046832 0.496284301 0.48130042 0.465318215 0.448052021 0.429086872 0.418766848 0.407731158 0.39581141 0.382745995 0.37392045 0.37355961 0.368150017 0.35810363 0.35768265 0.351340887 0.33952037 0.33891897 0.33098049 0.32851475 0.32773293 0.31582521 0.31579514 0.31489304 0.30352658 0.300104614 0.29880559 0.296322409 0.28803632 0.288004446 0.286684373 0.278054884 0.276346307 0.262565226 0.252503804 0.235375932 0.272034 0.264902 0.257857 0.251177 0.245477 0.241946 0.237573 0.231703 0.227394 0.224471];    

DIPPR_exp = DIPPR_exp(Temp_exp<=275);
Temp_exp = Temp_exp(Temp_exp<=275);
DIPPR_exp = DIPPR_exp(Temp_exp>=175);
Temp_exp = Temp_exp(Temp_exp>=175);

rhoc_g = 0.2015;
TC_g = 305.347;
A_g = 0.00048;
B_g = 0.05697;
beta_g = 0.33613;

[rhoc_fit, TC_fit, A_fit, B_fit, beta_fit] = fit_rho_L_beta(Temp_exp,DIPPR_exp,rhoc_g, TC_g, A_g, B_g, beta_g);

n = length(Temp_exp);

p = 5;

rho_L_hat = @(T,rhoc,TC,A,B,beta) rhoc + A .* (TC - T) + B .* (TC - T) .^ beta;

SE_fit = ((DIPPR_exp - rho_L_hat(Temp_exp,rhoc_fit,TC_fit,A_fit,B_fit,beta_fit))).^2;
SSE_fit = sum(SE_fit);

sigma = SSE_fit/(n-p);

RHS = sigma * (n + p * (finv(0.95,p,n-p)-1)); % No 0.95^p because I want the true 95%

% When using the entire temperature range
% A_range = linspace(0.935,1.065,20)*A_fit;
% B_range = linspace(0.95,1.06,20)*B_fit;
% rhoc_range = linspace(0.975,1.025,10)*rhoc_fit;
% TC_range = linspace(0.9999,1.00015,10)*TC_fit;
% beta_range = linspace(0.965,1.035,40)*beta_fit;

% When using the range of simulations

A_range = linspace(0.95,1.05,20)*A_fit;
B_range = linspace(0.945,1.06,20)*B_fit;
rhoc_range = linspace(0.96,1.04,20)*rhoc_fit;
TC_range = linspace(0.997,1.003,20)*TC_fit;
beta_range = linspace(0.97,1.03,40)*beta_fit;

Ext = zeros(1,(length(A_range)*length(B_range)*length(rhoc_range)*length(TC_range)*length(beta_range)));
A_ext = Ext;
B_ext = Ext;
rhoc_ext = Ext;
TC_ext = Ext;
beta_ext = Ext;
clear Ext

s=1;

for g=1:length(A_range)
    
    for h=1:length(B_range)

        for i=1:length(rhoc_range)
    
            for j=1:length(TC_range)
                
                for k=1:length(beta_range)
      
       SE = ((DIPPR_exp - rho_L_hat(Temp_exp,rhoc_range(i),TC_range(j),A_range(g),B_range(h),beta_range(k)))).^2;
        
       SSE = sum(SE);
       
       if SSE < RHS 
                 
           A_ext(s) = A_range(g);
           B_ext(s) = B_range(h);
           rhoc_ext(s) = rhoc_range(i);
           TC_ext(s) = TC_range(j);
           beta_ext(s) = beta_range(k);
           
           s=s+1;
           
           if SSE < SSE_fit % This is to make sure that we actually have the global minimum. If not, rerun the analysis with the new optimum (after plugging in as a new Mathcad guess)
               
              new_best_fit = [A_range(g) B_range(h) rhoc_range(i) TC_range(j) beta_range(k)]; 
              SSE_fit = SSE;
              
           end
           
       end
       
                end
       
            end
            
       end
       
   end
    
end

% This eliminates the superfluous zero elements
A_ext = A_ext(:,1:(s-1));
B_ext = B_ext(:,1:(s-1));
rhoc_ext = rhoc_ext(:,1:(s-1));
TC_ext = TC_ext(:,1:(s-1));
beta_ext = beta_ext(:,1:(s-1));

A_low = min(A_ext);
A_high = max(A_ext);
B_low = min(B_ext);
B_high = max(B_ext);
rhoc_low = min(rhoc_ext);
rhoc_high = max(rhoc_ext);
TC_low = min(TC_ext);
TC_high = max(TC_ext);
beta_low = min(beta_ext);
beta_high = max(beta_ext);

A_low/A_fit
A_high/A_fit
B_low/B_fit
B_high/B_fit
rhoc_low/rhoc_fit
rhoc_high/rhoc_fit
TC_low/TC_fit
TC_high/TC_fit
beta_low/beta_fit
beta_high/beta_fit

subplot(5,2,1)
plot(A_ext,B_ext)
subplot(5,2,2)
plot(rhoc_ext,B_ext)
subplot(5,2,3)
plot(A_ext,TC_ext)
subplot(5,2,4)
plot(rhoc_ext,TC_ext)
subplot(5,2,5)
plot(A_ext,rhoc_ext)
subplot(5,2,6)
plot(B_ext,TC_ext)
subplot(5,2,7)
plot(A_ext,beta_ext)
subplot(5,2,8)
plot(B_ext,beta_ext)
subplot(5,2,9)
plot(rhoc_ext,beta_ext)
subplot(5,2,10)
plot(TC_ext,beta_ext)


T_scan = linspace(min(Temp_exp), max(Temp_exp),1000);

for i = 1:length(T_scan)
    
   rho_L = rho_L_hat(T_scan(i), rhoc_ext, TC_ext, A_ext, B_ext, beta_ext);
   rho_L_max(i) = max(rho_L);
   rho_L_min(i) = min(rho_L);

end

rho_L = rho_L_hat(T_scan,rhoc_fit,TC_fit,A_fit,B_fit,beta_fit);

error_rho_L = (rho_L_max - rho_L_min)./rho_L;

figure
hold
scatter(Temp_exp,DIPPR_exp,'r')
plot(T_scan,rho_L)
plot(T_scan,rho_L_max,'--')
plot(T_scan,rho_L_min,'--')
hold

figure
plot(T_scan,error_rho_L)


