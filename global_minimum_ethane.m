% This shows the relationship between sigma, epsilon and rho_L

% If removing 45.2 from the octane runs
% liq_ave_ave = liq_ave_ave(:,2:5,:);
% vap_ave_ave = vap_ave_ave(:,2:5,:);
% epsilon = epsilon(2:5);

% Use this to just rerun this script (some variables might need to be
% deleted that haven't been if parameter space dimensions has changed)
% clear liq_hat liq_dev liq_acceptable liq_opt acceptable Y_sigma_refined Y_global RMS_opt TC_opt rhoc_opt liq_hat_plot liq_acceptable_plot SSE_opt

% Use this if you don't intend on revisiting the previous script
clearvars -except liq_ave_ave vap_ave_ave n_sim Temp_sim Temp_exp DIPPR_exp epsilon sigma DIPPR_TC MW LDN_Funke LDN_Funke_sim LDN_DIPPR_exp LDN_DIPPR_sim

sigma_refined = linspace(min(sigma),max(sigma),50);
epsilon_refined = linspace(min(epsilon),max(epsilon),50);

% For plotting comparison of ethane
% Type B
% sigma_refined = linspace(3.72,3.78,50);
% epsilon_refined = linspace(97,100,50);

% For plotting comparison of ethane
% Type A
% sigma_refined = linspace(3.7471,3.751,20);
% epsilon_refined = linspace(98.38,98.62,20);

% With the more refined simulation data
% sigma_refined = linspace(3.746,3.7515,500);
% epsilon_refined = linspace(98.32,98.65,500);

% With the more refined simulation data but at the 99% confidence
% sigma_refined = linspace(3.746,3.752,100);
% epsilon_refined = linspace(98.32,98.68,100);

% With the more refined simulation data but with Type B analysis
% sigma_refined = linspace(3.73,3.77,200);
% epsilon_refined = linspace(97.5,99.5,200);

% With the more refined simulation data but with Type B analysis
% sigma_refined = linspace(3.72,3.78,200);
% epsilon_refined = linspace(97.3,99.7,200);

% With the more refined simulation data but with true Type B analysis
% Assigning entire correlation uncertainty
% sigma_refined = linspace(3.72,3.78,100);
% epsilon_refined = linspace(97.3,100.3,100);

% With the more refined simulation data but with true Type B analysis
% Assigning entire correlation uncertainty
% Using the DIPPR correlation instead
% sigma_refined = linspace(3.71,3.78,100);
% epsilon_refined = linspace(97.1,100,100);

% Doing octane (I should just make it's own file)
% sigma_refined = linspace(3.9,3.985,20);
% epsilon_refined = linspace(45,45.9,20);

% Doing octane with the more refined short simulations and just TRC data
% sigma_refined = linspace(3.96,3.99,50);
% epsilon_refined = linspace(45.25,45.7,50);

% Doing octane with the more refined short simulations and just the data I accepted
% sigma_refined = linspace(3.97,3.983,30);
% epsilon_refined = linspace(45.35,45.6,30);

% Doing octane with the more refined short simulations and just the data I
% accepted but at 99%
% sigma_refined = linspace(3.967,3.986,400);
% epsilon_refined = linspace(45.33,45.63,400);

% For TraPPE's CH3
% Doing octane with the more refined short simulations but with 8 replicates 
% and just the data I accepted but at 99%
% sigma_refined = linspace(3.966,3.981,30);
% epsilon_refined = linspace(45.38,45.60,30);

% Doing octane with the more refined short simulations but with 8 replicates 
% and just the data I accepted but at 99%
% sigma_refined = linspace(3.966,3.981,400);
% epsilon_refined = linspace(45.38,45.60,400);

% For my CH3
% Doing octane with the more refined short simulations but with 8 replicates 
% and just the data I accepted but at 99%
% sigma_refined = linspace(3.967,3.982,40);
% epsilon_refined = linspace(45.26,45.50,40);

% For non-optimal CH3, just the same region as before for graphical
% purposes
% sigma_refined = linspace(3.966,3.983,30);
% epsilon_refined = linspace(45.20,45.55,30);

% In order to compare all the CH3 ones, I need to recenter this using the
% same range of CH2 simulations
% sigma_refined = linspace(3.961,3.982,30);
% epsilon_refined = linspace(45.25,45.5,30);

% Doing octane with the more refined short simulations but with 8 replicates 
% without the highest temperature value
% sigma_refined = linspace(3.965,3.979,40);
% epsilon_refined = linspace(45.3,45.46,40);

% Doing octane with the more refined short simulations but with 8 replicates 
% without the highest temperature value at 99%
% sigma_refined = linspace(3.963,3.981,100);
% epsilon_refined = linspace(45.28,45.48,100);

% Doing octane with the more refined short simulations but with 8 replicates 
% without the highest temperature value at 99% for a Type B analysis
% sigma_refined = linspace(3.935,4.01,300);
% epsilon_refined = linspace(44.95,45.8,300);

% Doing octane with the more refined short simulations but with 8 replicates 
% without the highest temperature value at 99% for a true Type B analysis
% The uncertainty is for the entire correlation
% sigma_refined = linspace(3.92,4.02,300);
% epsilon_refined = linspace(44.6,46.1,300);

% Doing octane with the more refined short simulations but with 8 replicates 
% without the highest temperature value at 99% for a true Type B analysis
% The uncertainty is for the entire correlation
% Using the DIPPR correlation instead
% sigma_refined = linspace(3.9,4.02,100);
% epsilon_refined = linspace(44.5,46.1,100);

for g = 1:2

for j = 1:length(Temp_sim)
    
    X = ones(length(sigma),1);
    X = [X, sigma'];
    X_refined = ones(length(sigma_refined),1);
    X_refined = [X_refined, sigma_refined'];
    
% figure

for h = 1:length(epsilon) % Span the different epsilons
    
    if g == 1
        
        Y = vap_ave_ave(j,h,:);
        
    else
        
        Y = liq_ave_ave(j,h,:);
        
    end
    
    Y = permute(Y,[3 2 1]);
    
    b_fit = (X'*X)\(X'*Y);
    
    Y_fit = X*b_fit;
    Y_refined = X_refined*b_fit;
    
    Y_sigma(j,h,:) = Y_fit;
    Y_sigma_refined(j,h,:) = Y_refined;
    
%     hold
%     scatter(sigma,Y)
%     plot(sigma,Y_fit)
%     title(Temp_sim(j))
%     hold
    
end

%     for h = 1:length(epsilon_refined)
%    
%         Y_sigma_refined(j,h,:) = 0;
%     
%     end

    X = ones(length(epsilon),1);
    X = [X, epsilon'];
    X_refined = ones(length(epsilon_refined),1);
    X_refined = [X_refined, epsilon_refined'];
    
    for i = 1:length(sigma_refined)
        
    Y_sig_refined = Y_sigma_refined(j,:,i);
    Y_sig_refined = permute(Y_sig_refined,[2 3 1]);
    
    b_fit_refined = (X'*X)\(X'*Y_sig_refined);
    Y_fit_refined = X_refined*b_fit_refined;
    
    Y_global(j,:,i) = Y_fit_refined;
    
%     hold
%     scatter(epsilon,Y_sig_refined,'g')
%     plot(epsilon_refined,Y_fit_refined,'r')
%     title(Temp_sim(j))
%     hold
    
    end

% figure

for i = 1:length(sigma)
    
    if g == 1
        
        Y = vap_ave_ave(j,:,i);
        
    else
        
        Y = liq_ave_ave(j,:,i);
        
    end
    
    Y = permute(Y,[2 3 1]);
    
    Y_sig = Y_sigma(j,:,i);
    Y_sig = permute(Y_sig,[2 3 1]);
        
    b_fit = (X'*X)\(X'*Y);
    b_eps_sig = (X'*X)\(X'*Y_sig);
    
    Y_fit = X*b_fit;
    Y_eps_sig = X*b_eps_sig;
    
    Y_epsilon(j,:,i) = Y_fit;
    Y_epsilon_sigma(j,:,i) = Y_eps_sig;
    
%     hold
%     scatter(epsilon,Y,'b')
%     scatter(epsilon,Y_sig,'g')
%     plot(epsilon,Y_fit,'b')
%     plot(epsilon,Y_eps_sig,'g')
%     title(Temp_sim(j))
%     hold
    
end

end

    if g == 1

        vap_smooth = Y_global;
    
    else
        
        liq_smooth = Y_global;
        
    end

end

if MW == 30.07

    beta = 0.326; % Seeing what the beta effect is
   
else
    
    beta = 0.32;
    
end
  
   bb = 1; % Not used here unless modify code for different betas
    
    TC = zeros(length(epsilon_refined),length(sigma_refined));
    rhoc = TC;
    rhoc_fit = TC;
    TC_fit = TC;
    PC = TC;
    
if MW == 226.45 % If doing C16 you don't need all this extra stuff
    
   for h = 1:length(epsilon_refined) % Cycles through all the different epsilons
       
        for i = 1:length(sigma_refined) % Cycles through all the different sigmas 
            
            [TC_opt(h,i), rhoc_opt(h,i), A, b] = towhee_error_model(Temp_sim,vap_smooth(:,h,i)',liq_smooth(:,h,i)',beta); % This uses my error model
            
        end
        
   end
   
else
    
       Temp_plot = linspace(min(Temp_exp),max(Temp_exp),1000);
%    Temp_plot = linspace(0.5*min(min(TC_fit)),min(min(TC_fit)),1000); % Have to do this after you know what TC_fit is from doing it first with the above Temp_plot
    
valid95 = zeros(length(epsilon_refined),length(sigma_refined));

U_DIPPR = 0.01;
U_abs_DIPPR = U_DIPPR*mean(LDN_DIPPR_exp); % Using the average of experimental data
% U_abs_DIPPR = U_DIPPR*max(LDN_DIPPR_exp); % Using the maximum absolute error

for h = 1:length(epsilon_refined) % Cycles through all the different epsilons
       
   for i = 1:length(sigma_refined) % Cycles through all the different sigmas
       
   [TC(h,i), rhoc(h,i), A, b] = towhee_error_model(Temp_sim,vap_smooth(:,h,i)',liq_smooth(:,h,i)',beta); % This uses my error model
     
   liq_TR = rhoc(h,i) + A*TC(h,i)*(0.3) + b*TC(h,i)^(beta)*(0.3)^beta;
   
   [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L(Temp_sim,liq_smooth(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta); % This fits the smoothed rho_L data with the weighting model
     
%    [rhoc_fit(h,i), TC_fit(h,i), A, b] = fit_rho_L_constrained(Temp,liq_smooth(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta,max(Temp_exp)); % This just fits the rho_L data but keeps TC above the maximum experimental temperature without weighting the data
    
% % % Using Rackett equation

%      PC_g = PC_from_Rackett_alt(TC(h,i),rhoc(h,i),liq_TR,0.7,MW);
%      
%      [rhoc_fit(h,i), ZC_fit(h,i), TC_fit(h,i), D_fit(h,i)] = weighted_fit_rho_L_Rackett(Temp_sim,liq_smooth(:,h,i)',rhoc(h,i),TC(h,i),PC_g,MW);
% 
%      rho_L_Rackett = @(T,b) b(1) .* b(2) .^ (-(1-T./b(3)).^(b(4)));
%      
%      rho_L_Rackett_fit = @(T) rho_L_Rackett(T,[rhoc_fit(h,i), ZC_fit(h,i), TC_fit(h,i), D_fit(h,i)]);
%      
%      liq_hat(:,h,i) = rho_L_Rackett_fit(Temp_exp);
%      
%      liq_hat_sim(:,h,i) = rho_L_Rackett_fit(Temp_exp);
%      
%      liq_hat_plot(:,h,i) = rho_L_Rackett_fit(Temp_plot);
%           
%      liq_dev(:,h,i) = rho_L_Rackett_fit(Temp_sim); - liq_smooth(:,h,i)';

% % % End of Rackett section
     
     liq_hat(:,h,i) = rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_exp) + b.*(TC_fit(h,i)-Temp_exp).^(beta);
     
     liq_hat_sim(:,h,i) = rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_sim) + b.*(TC_fit(h,i)-Temp_sim).^(beta);
      
     PC(h,i) = PC_from_Rackett_alt(TC(h,i),rhoc(h,i),liq_hat_sim(1,h,i),Temp_sim(1)/TC(h,i),MW);
        
%      liq_hat_plot(:,h,i) = rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_plot) + b.*(TC_fit(h,i)-Temp_plot).^(beta);
%           
%      liq_dev(:,h,i) = (rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_sim) + b.*(TC_fit(h,i)-Temp_sim).^(beta)) - liq_smooth(:,h,i)';
% 
%      SSE_dev(h,i) = sum(liq_dev(:,h,i).^2);
%      
%      RMS_dev(h,i) = sqrt(SSE_dev(h,i)/n_sim);
% 
%      per_dev(:,h,i) = liq_dev(:,h,i) ./ liq_smooth(:,h,i) * 100;
%      
%      error_hat = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(Temp_exp/DIPPR_TC));

     % This is the true Type B analysis on the correlation uncertaint
       
%      Dev_hi = DIPPR_exp'*(1+U_DIPPR) - liq_hat(:,h,i);  
%      Dev_lo = DIPPR_exp'*(1-U_DIPPR) - liq_hat(:,h,i);
     
%      Dev_hi = LDN_Funke'*(1+U_DIPPR) - liq_hat(:,h,i);  
%      Dev_lo = LDN_Funke'*(1-U_DIPPR) - liq_hat(:,h,i);
     
     Dev_hi = LDN_DIPPR_exp'*(1+U_DIPPR) - liq_hat(:,h,i);  
     Dev_lo = LDN_DIPPR_exp'*(1-U_DIPPR) - liq_hat(:,h,i);

    % Dr. Giles thinks a constant absolute error might be better
    
%      Dev_hi = LDN_DIPPR_exp' + U_abs_DIPPR - liq_hat(:,h,i);  
%      Dev_lo = LDN_DIPPR_exp' - U_abs_DIPPR - liq_hat(:,h,i);
% %      
     % I really could just look at the temperatures simulated since I am
     % using a DIPPR correlation here
     
%      Dev_hi = LDN_DIPPR_sim'*(1+U_DIPPR) - liq_hat_sim(:,h,i);  
%      Dev_lo = LDN_DIPPR_sim'*(1-U_DIPPR) - liq_hat_sim(:,h,i);

    % Dr. Giles thinks a constant absolute error might be better
    
%      Dev_hi = LDN_DIPPR_sim' + U_abs_DIPPR - liq_hat_sim(:,h,i);  
%      Dev_lo = LDN_DIPPR_sim' - U_abs_DIPPR - liq_hat_sim(:,h,i);
      
     
     if min(Dev_hi) >= 0 && max(Dev_lo) <= 0
        
         valid95(h,i) = 1;
         
     end
     
     if TC_fit(h,i) > max(Temp_exp)
   
        SSE(h,i) = sum((DIPPR_exp'-liq_hat(:,h,i)).^2); % Non-weighted approach
%         SSE(h,i) = sum((DIPPR_exp'-liq_ave(:,h,i)).^2); % Non-weighted approach, using the average of the replicates (only valid when using TRC correlation)
%         SSE(h,i) = sum(((DIPPR_exp'-liq_hat(:,h,i))./error_hat').^2); % Weighted approach

%         SSE(h,i) = SSE(h,i) + SSE_dev(h,i);
         
     else
         
        SSE(h,i) = 100000000;
         
     end
     
%         SSE(h,i) = sum((TRC'-liq_ave(:,h,i)).^2); % Non-weighted approach, using all the replicate values
%         SSE(h,i) = sum(((TRC'-liq_ave(:,h,i))./error_hat').^2); % Weighted approach, using all the replicate values

        
     if bb == 1
               
                SSE_opt(h,i) = SSE(h,i);
                beta_opt(h,i) = beta;
                TC_opt(h,i) = TC(h,i);
                rhoc_opt(h,i) = rhoc(h,i);
                liq_opt(:,h,i) = liq_hat(:,h,i);
                
     elseif SSE(h,i) < SSE_opt(h,i)
            
                SSE_opt(h,i) = SSE(h,i);
                beta_opt(h,i) = beta;
                TC_opt(h,i) = TC(h,i);
                rhoc_opt(h,i) = rhoc(h,i);
                liq_opt(:,h,i) = liq_hat(:,h,i);

     end
     
   end
   
end

    VC = MW ./ rhoc_opt;

    ZC = PC .* VC ./ TC / 8.314472;

% Trying to automate the contours
    p = 2; % Not sure what this value really should be, depends on the method
    n = length(Temp_exp);

    alpha = 0.95;

    [SSE_min, I] = min(SSE_opt);
    [SSE_min, J] = min(SSE_min);
    
    eps_opt = epsilon_refined(I(J));
    sig_opt = sigma_refined(J);
    opt = [eps_opt, sig_opt];
    
    RHS = SSE_min * (1 + (p/(n-p))*finv(alpha,p,n-p));
    
    sigma2 = SSE_min/(n-p);
    
    sigma2 = 0.002403^2; % This is from the DIPPR uncertainty for ethane
    sigma2 = 0.00278^2; % This is from the DIPPR uncertainty for octane
    
    PDF = fcdf((SSE_opt - SSE_min)/(sigma2*p),p,n-p);
    
    RMS_opt = sqrt(SSE_opt/n);
%     RMS_min = sqrt(SSE_min/n);
%     RHS_RMS = sqrt(RHS/n);
%     
    z = 1;
    
    figure
    hold
    scatter(Temp_exp,DIPPR_exp,'r')
%     plot(Temp_plot',liq_hat_plot(:,I(J),J)','b')
    hold
           
for h = 1:length(epsilon_refined)
    
    for i = 1:length(sigma_refined)

        if valid95(h,i) == 1 
%         if SSE_opt(h,i) <= RHS
%       Pareto Front analysis for ethane
%         if valid95(h,i) == 1 && TC_opt(h,i) < 305.93 && TC_opt(h,i) > 304.7 % Using just experimental uncertainty
%         if valid95(h,i) == 1 && TC_opt(h,i) < (305.93+0.8) && TC_opt(h,i) > (304.7-0.8) % Using experimental, numerical, and surrogate model uncertainty
%       Pareto Front analysis for octane
%         if valid95(h,i) == 1 && TC_opt(h,i) < 569.84 && TC_opt(h,i) > 567.56 % Using just experimental uncertainty
%         if valid95(h,i) == 1 && TC_opt(h,i) < (569.84+2.8) && TC_opt(h,i) > (567.56-2.8) % Using experimental, numerical, and surrogate model uncertainty
%             hold
%             plot(Temp_plot',liq_hat_plot(:,h,i)','g')
% %             scatter(Temp,liq_ave(:,h,i)','g') % This doesn't work because we need h and i to span just the simulations to use this 
%             hold
            acceptable(:,z) = [epsilon_refined(h) sigma_refined(i)];
            TC_acceptable(z) = TC_opt(h,i);
            rhoc_acceptable(z) = rhoc_opt(h,i);
            PC_acceptable(z) = PC(h,i);
            liq_acceptable(:,z) = liq_hat(:,h,i);
            liq_acceptable_sim(:,z) = liq_hat_sim(:,h,i);
            z=z+1;
%             hold
% % %             scatter(Temp_sim,liq_dev(1:n_sim,h,i)','r')
% %             scatter(Temp_sim,per_dev(1:n_sim,h,i)','r')
% % %             plot(Temp_sim,std_sim(:,h,i))
% % %             plot(Temp_sim,-std_sim(:,h,i))
% %             plot(Temp_sim,per_error_95_sim(:,h,i),'g')
% %             plot(Temp_sim,-per_error_95_sim(:,h,i),'g')
%             xlabel('Temperature (K)')
%             ylabel('Percent Deviation')
%             hold
        end
        
%         scatter(Temp_sim,std_sim(:,h,i))
    
    end
    
end

[liq_hi, z_max] = max(liq_acceptable');
[liq_lo, z_min] = min(liq_acceptable');
liq_opt = liq_hat(:,I(J),J);

[liq_hi_sim, z_max] = max(liq_acceptable_sim');
[liq_lo_sim, z_min] = min(liq_acceptable_sim');

hold
plot(Temp_exp,liq_hi,'g')
plot(Temp_exp,liq_lo,'g')
plot(Temp_sim,liq_hi_sim,'c')
plot(Temp_sim,liq_lo_sim,'c')
% plot(Temp_sim,LDN_DIPPR_sim,'b')
% plot(Temp_sim,LDN_DIPPR_sim*(1+U_DIPPR),'b--')
% plot(Temp_sim,LDN_DIPPR_sim*(1-U_DIPPR),'b--')
hold
% 
% extrema_high = acceptable(:,z_max(:));
% extrema_low = acceptable(:,z_min(:));
% 
% dev_hi = (liq_hi - DIPPR_exp)./DIPPR_exp * 100;
% dev_lo = (liq_lo - DIPPR_exp)./DIPPR_exp * 100;
% 
% [liq_hi_plot, z_max_plot] = max(liq_acceptable_plot');
% [liq_lo_plot, z_min_plot] = min(liq_acceptable_plot');
% liq_opt_plot = liq_hat_plot(:,I(J),J);
% 
% extrema_high_plot = acceptable(:,z_max_plot(:));
% extrema_low_plot = acceptable(:,z_min_plot(:));
% 
% % figure
% % hold
% % scatter(extrema_high_plot(2,:),extrema_high_plot(1,:),'b')
% % scatter(extrema_low_plot(2,:),extrema_low_plot(1,:),'g')
% % scatter(extrema_high(2,:),extrema_high(1,:),'rx')
% % scatter(extrema_low(2,:),extrema_low(1,:),'kx')
% % hold
% 
% liq_width = liq_hi_plot-liq_lo_plot;
% 
% per_uncertainty = (liq_width')./ liq_opt_plot * 100;
% 
% % figure
% % hold
% % scatter(Temp_exp,dev_hi)
% % scatter(Temp_exp,dev_lo)
% % plot(Temp_plot,per_uncertainty)
% % plot(Temp_plot,-per_uncertainty)
% % hold
% 
% eps_max = max(max(acceptable'));
% sig_max = min(max(acceptable'));
% eps_min = max(min(acceptable'));
% sig_min = min(min(acceptable'));
% 
% eps_uncertainty = (eps_max - eps_min) / eps_opt / 2 * 100;
% sig_uncertainty = (sig_max - sig_min) / sig_opt / 2 * 100;
% 
% eps_STD = (eps_max - eps_min) / 2 / 1.96 * sqrt(n);
% sig_STD = (sig_max - sig_min) / 2 / 1.96 * sqrt(n);
% 
% eps_per_STD = eps_STD / eps_opt * 100;
% sig_per_STD = sig_STD / sig_opt * 100;
% 
figure
contour(sigma_refined,epsilon_refined,RMS_opt)
% 
% [TC_hi, z_TC_max] = max(TC_acceptable);
% extrema_high_TC = acceptable(:,z_TC_max);
% [TC_lo, z_TC_min] = min(TC_acceptable);
% extrema_low_TC = acceptable(:,z_TC_min);
% 
% [rhoc_hi, z_rhoc_max] = max(rhoc_acceptable);
% extrema_high_rhoc = acceptable(:,z_rhoc_max);
% [rhoc_lo, z_rhoc_min] = min(rhoc_acceptable);
% extrema_low_rhoc = acceptable(:,z_rhoc_min);
% 
% figure
% contour(sigma_refined,epsilon_refined,TC_opt)
% 
% figure
% contour(sigma_refined,epsilon_refined,rhoc_opt)
% 
figure
hold
scatter(acceptable(2,:),acceptable(1,:),'b')
scatter(sig_opt,eps_opt,'r')
hold
% 
% % This way isn't as rigorous, but it appeared to work just as well
% 
% % for i = 1:10000
% %     % r=normcdf(normrnd(0,1),0,1)*ones(length(epsilon_refined),length(sigma_refined));
% %     r=fcdf(frnd(p,n-p),p,n-p)*ones(length(epsilon_refined),length(sigma_refined));
% %     dif = abs(r-PDF);
% %     [dif_min, I] = min(dif);
% %     [dif_min, J] = min(dif_min);
% %     eps_norm(i) = epsilon_refined(I(J));
% %     sig_norm(i) = sigma_refined(J); 
% %     TC_norm(i) = TC_opt(I(J),J);
% %     rhoc_norm(i) = rhoc_opt(I(J),J);
% % end
% 
% histogram_sampling
% 
% [TC_counts, TC_centers] = hist(TC_norm,1000);
% 
% [TC_lo_rig, TC_hi_rig, actual_confidence_TC] = integrate_histogram(TC_counts,TC_centers,alpha);
% 
% TC_uncertainty = (TC_hi_rig - TC_lo_rig) / (TC_hi_rig + TC_lo_rig) * 100;
% 
% [rhoc_counts, rhoc_centers] = hist(rhoc_norm,1000);
% 
% [rhoc_lo_rig, rhoc_hi_rig, actual_confidence_rhoc] = integrate_histogram(rhoc_counts,rhoc_centers,alpha);
% 
% rhoc_uncertainty = (rhoc_hi_rig - rhoc_lo_rig)/(rhoc_hi_rig + rhoc_lo_rig)*100;
% 
% [PC_counts, PC_centers] = hist(PC_norm,1000);
% 
% [PC_lo_rig, PC_hi_rig, actual_confidence_PC] = integrate_histogram(PC_counts,PC_centers,alpha);
% 
% PC_uncertainty = (PC_hi_rig - PC_lo_rig) / (PC_hi_rig + PC_lo_rig) * 100;
% 
% [ZC_counts, ZC_centers] = hist(ZC_norm,1000);
% 
% [ZC_lo_rig, ZC_hi_rig, actual_confidence_ZC] = integrate_histogram(ZC_counts,ZC_centers,alpha);
% 
% ZC_uncertainty = (ZC_hi_rig - ZC_lo_rig) / (ZC_hi_rig + ZC_lo_rig) * 100;
% % % % 
% % % % for k = 1:length(Temp_exp)
% % % %    
% % % %     [rhoL_counts(:,k), rhoL_centers(:,k)] = hist(rhoL_norm(:,k),1000);
% % % %         
% % % %     [rhoL_lo_rig(k), rhoL_hi_rig(k)] = integrate_histogram(rhoL_counts(:,k),rhoL_centers(:,k),alpha);
% % % %     
% % % %     rhoL_uncertainty(k) = (rhoL_hi_rig(k) - rhoL_lo_rig(k)) / (rhoL_hi_rig(k) + rhoL_lo_rig(k)) * 100;
% % % %     
% % % % end
% % 
% [eps_counts, eps_centers] = hist(eps_norm,1000);
% 
% [eps_lo_rig, eps_hi_rig, actual_confidence_eps] = integrate_histogram(eps_counts,eps_centers,alpha);
% 
% eps_uncertainty = (eps_hi_rig - eps_lo_rig) / (eps_hi_rig + eps_lo_rig) * 100;
% 
% [sig_counts, sig_centers] = hist(sig_norm,1000);
% 
% [sig_lo_rig, sig_hi_rig, actual_confidence_sig] = integrate_histogram(sig_counts,sig_centers,alpha);
% 
% sig_uncertainty = (sig_hi_rig - sig_lo_rig) / (sig_hi_rig + sig_lo_rig) * 100;
% 
% figure
% subplot(2,2,1)
% hist(TC_norm,1000)
% subplot(2,2,2)
% hist(rhoc_norm,1000)
% subplot(2,2,3)
% hist(eps_norm,1000)
% subplot(2,2,4)
% hist(sig_norm,1000)
% 
% figure
% hist(PC_norm,1000)
% % 
% figure
%     for k = 1:length(Temp_exp)
%         subplot(3,4,k)
%         hist(rhoL_norm(:,k),1000)
%     end

liq_ave_ave_contour = permute(liq_ave_ave,[2 3 1]);
liq_smooth_contour = permute(liq_smooth,[2 3 1]);
liq_hat_contour = permute(liq_hat_sim,[2 3 1]);

for j = 1:length(Temp_sim)
    
figure
subplot(1,3,1)
contour(sigma, epsilon, liq_ave_ave_contour(:,:,j))
title(['Raw data'])
xlabel('Sigma (Ang)')
ylabel('Epsilon (K)')
subplot(1,3,2)
contour(sigma_refined, epsilon_refined, liq_smooth_contour(:,:,j))
title(['Smoothed'])
xlabel('Sigma (Ang)')
ylabel('Epsilon (K)')
subplot(1,3,3)
contour(sigma_refined, epsilon_refined, liq_hat_contour(:,:,j))
title(['After Regression to LDN Model'])
xlabel('Sigma (Ang)')
ylabel('Epsilon (K)')

end
% % 
% % border_eps_sig

end