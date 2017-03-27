Mw_ethane = 30.07;

b_LJ_fit = [0.21729442062711482 0.05693608598247599 3.091216622550399 0.2703919613194321 0.31936400699635553]; % This is not weighting the simulation results
b_LJ_fit = [0.2173805215340362 0.05698297259417123 3.0910955219325817 0.270275670063868 0.31929819333035586]; % This is weighted
b_LJ_fit = [0.21730244863840334 0.05705205696124398 3.0911267468395316 0.27029597123869437 0.31913968475996]; % This is using the more refined grid and a weighted optimization
b_LJ_fit = [0.21724452840211333 0.05704063825248595 3.091275216548376 0.27029104539701804 0.31925075168064304]; % This is using only the range of simulations that were acceptable

rho_plus_LJ = @(T) b_LJ_fit(1) + b_LJ_fit(2).*(b_LJ_fit(3)-T) + b_LJ_fit(4).*(b_LJ_fit(3)-T).^b_LJ_fit(5);

rho_LJ = @(T,epsilon,sigma) rho_plus_LJ(T./epsilon) * Mw_ethane ./ (sigma.^3) * 0.0016605778811026233; % Converts this to gm/mL

% Optimal using nonweighted
eps_LJ_opt = 98.48300694436236;
sig_LJ_opt = 0.37491573562839176;

% Optimal using weighted
eps_LJ_opt = 98.48185817341388;
sig_LJ_opt = 0.3749148156724503;

% Optimal using the more refined simulations and weighted optimization
eps_LJ_opt = 98.48519716419385;
sig_LJ_opt = 0.37492074823253246;

% Optimal using only the range that encompassed the acceptable most refined in fitting and weighted
eps_LJ_opt = 98.48821093475672;
sig_LJ_opt = 0.374910620368508;

SSE_LJ = @(epsilon,sigma) sum((DIPPR_exp - rho_LJ(Temp_exp,epsilon,sigma)).^2);

sigma_refined = linspace(3.7471,3.751,100);
epsilon_refined = linspace(98.33,98.64,100);

SSE_global_min = SSE_LJ(eps_LJ_opt,sig_LJ_opt);

p = 2;

RHS_global = SSE_global_min * (1 + (p/(n-p))*finv(alpha,p,n-p));

SSE_global = zeros(length(epsilon_refined),length(sigma_refined));

z = 1;

T_plot = linspace(min(Temp_exp),max(Temp_exp),100);

sig_hi = zeros(length(epsilon_refined),1);
sig_lo = 1000*ones(length(epsilon_refined),1);
eps_acceptable = zeros(length(epsilon_refined),1);

for h = 1:length(epsilon_refined)
    
    for i = 1:length(sigma_refined)
        
        SSE_global(h,i) = SSE_LJ(epsilon_refined(h),sigma_refined(i)/10);
        
        if SSE_global(h,i) <= RHS_global
            
            acceptable(:,z) = [epsilon_refined(h), sigma_refined(i)];
            liq_acceptable(:,z) = rho_LJ(T_plot,epsilon_refined(h),sigma_refined(i)/10);
            liq_acceptable_exp(:,z) = rho_LJ(Temp_exp,epsilon_refined(h),sigma_refined(i)/10);
            z = z+1;
            
            if sigma_refined(i) > sig_hi(h)
                
                sig_hi(h) = sigma_refined(i);
                
            end
                
            if sigma_refined(i) < sig_lo(h)
               
                sig_lo(h) = sigma_refined(i);
                
            end
            
            eps_acceptable(h) = epsilon_refined(h);
            
            
        end

    end
    
end

sig_lo = sig_lo(eps_acceptable~=0);
sig_hi = sig_hi(eps_acceptable~=0);
eps_acceptable = eps_acceptable(eps_acceptable~=0);

RMS_global = sqrt(SSE_global/n);
RMS_global_min = sqrt(SSE_global_min/n);

RHS_RMS_global = sqrt(RHS_global/n);

liq_hi = max(liq_acceptable');
liq_lo = min(liq_acceptable');
liq_opt = rho_LJ(T_plot,eps_LJ_opt,sig_LJ_opt);
liq_hi_exp = max(liq_acceptable_exp');
liq_lo_exp = min(liq_acceptable_exp');

dev_hi = (liq_hi_exp - DIPPR_exp)./DIPPR_exp * 100;
dev_lo = (liq_lo_exp - DIPPR_exp)./DIPPR_exp * 100;

liq_width = liq_hi-liq_lo;

% figure
% plot(T_plot,liq_width)
% 
% figure
% hold
% plot(Temp_exp,dev_hi)
% plot(Temp_exp,dev_lo)
% hold
% 
% figure
% hold
% scatter(Temp_exp,DIPPR_exp)
% plot(T_plot,liq_hi,'--')
% plot(T_plot,liq_lo,'--')
% plot(T_plot,liq_opt)
% hold

% figure
hold
scatter(acceptable(2,:),acceptable(1,:),'b')
scatter(sig_LJ_opt*10,eps_LJ_opt,'r')
% plot(sig_hi,eps_acceptable,'g')
% plot(sig_lo,eps_acceptable,'g')
hold