clear

MW = 30.07;

Temp_exp = [88.71 90.37 91 94.26 95 99.82 100 105.37 110 110.93 116.48 120 122.04 127.59 130 133.15 140 150 160 170 175 185 195 200 210 220 230 240 250 260 265 270 275 280 283.2 283.2 285 288.19 288.19 290 293.18 293.18 295 295.67 295.67 298 298.17 298.17 300 300.66 300.66 301 302 302.16 302.16 303 303.16 304.15 304.65 305.15 303.5 304 304.4 304.7 304.9 305 305.1 305.2 305.25 305.275]; 
DIPPR_exp = [0.6519176 0.65152669 0.650847108 0.64611409 0.646468916 0.64040079 0.640987155 0.63450707 0.630014612 0.62882384 0.62290005 0.618981929 0.61700633 0.6110224 0.60786505 0.60509861 0.59661887 0.585207305 0.573594271 0.561722635 0.555675558 0.543319795 0.530558087 0.52399982 0.51046832 0.496284301 0.48130042 0.465318215 0.448052021 0.429086872 0.418766848 0.407731158 0.39581141 0.382745995 0.37392045 0.37355961 0.368150017 0.35810363 0.35768265 0.351340887 0.33952037 0.33891897 0.33098049 0.32851475 0.32773293 0.31582521 0.31579514 0.31489304 0.30352658 0.300104614 0.29880559 0.296322409 0.28803632 0.288004446 0.286684373 0.278054884 0.276346307 0.262565226 0.252503804 0.235375932 0.272034 0.264902 0.257857 0.251177 0.245477 0.241946 0.237573 0.231703 0.227394 0.224471];    
rho_L_max_min % Calls a script that has the data from the regression uncertainty analysis
TRC = [0.551945935 0.528130197 0.500847103 0.472019478 0.437027813 0.396053755]; % TRC correlation
TraPPE = [0.551 0.527 0.499 0.469 0.432 0.396];
TraPPE_validation = [0.5512 0.5262 0.4984 0.4342 0.3937 0.3835 0.3726 0.3589];
Temp_validation = [178 197 217 256 275 279 283 288];
Temp = [178 197 217 236 256 275];
Temp_sim = Temp;
n_sim = length(Temp_sim);
[TC_exp, I_exp] = max(Temp_exp);
rhoc_exp = DIPPR_exp(I_exp);
DIPPR_TC = 305;
dLDN = 0.0012;
dTC = 0.04;
beta = 0.326;

DIPPR_high = DIPPR_exp(Temp_exp>=300);
Temp_high = Temp_exp(Temp_exp>=300);

% DIPPR_exp = DIPPR_exp(Temp_exp<=300);
% Temp_exp = Temp_exp(Temp_exp<=300);
DIPPR_exp = DIPPR_exp(Temp_exp<=275);
Temp_exp = Temp_exp(Temp_exp<=275);
DIPPR_exp = DIPPR_exp(Temp_exp>=178);
Temp_exp = Temp_exp(Temp_exp>=178);
% DIPPR_exp = [DIPPR_exp, rhoc_exp];
% Temp_exp = [Temp_exp, TC_exp];
% DIPPR_exp = [DIPPR_exp, DIPPR_high];
% Temp_exp = [Temp_exp, Temp_high];

% DIPPR_exp = TRC;
% Temp_exp = Temp;

Parameters_A = [0.205807539 305.2997254 0.000557459 0.055573024 0.330088978]; % Fit to just Funke's rho_r and rho_s data
Parameters_B = [0.201491585 305.3473084 0.000476052 0.056966233 0.336131842]; % Fit to all the rho_L data
Parameters_C = [0.197622849 305.3463998 0.000414702 0.058572282 0.337994322]; % Fit to just the rho_L from 185-275

Funke_fit = Parameters_C;

LDN_Funke = Funke_fit(1) + Funke_fit(3) .* (Funke_fit(2) - Temp_exp) + Funke_fit(4) .* (Funke_fit(2) - Temp_exp) .^ Funke_fit(5);
LDN_Funke_Sim = Funke_fit(1) + Funke_fit(3) .* (Funke_fit(2) - Temp_sim) + Funke_fit(4) .* (Funke_fit(2) - Temp_sim) .^ Funke_fit(5);

SSE_Funke = sum((DIPPR_exp'-LDN_Funke').^2); % Non-weighted approach

 p = 5; % This is just for Funke, not LJ
 n = length(Temp_exp);
 alpha = 0.95;

 RHS_Funke = SSE_Funke * (1 + (p/(n-p))*finv(alpha,p,n-p));
 
%  Parameters_DIPPR = [1.9122 0.27937 305.32 0.29187]; % These are the parameters that DIPPR has
 Parameters_DIPPR = [1.92794 0.2816696 305.32 0.2858483]; % These are the parameters if you fix TC and VC to those of DIPPR 
 
 LDN_DIPPR_exp = MW / 1000 * Parameters_DIPPR(1) ./ (Parameters_DIPPR(2) .^ (1+(1-Temp_exp/Parameters_DIPPR(3)).^Parameters_DIPPR(4)));
 LDN_DIPPR_sim = MW / 1000 * Parameters_DIPPR(1) ./ (Parameters_DIPPR(2) .^ (1+(1-Temp_sim/Parameters_DIPPR(3)).^Parameters_DIPPR(4)));
  
hold
scatter(Temp_exp,DIPPR_exp)
scatter(Temp,TRC,'r')
% % plot(Temp,TRC,'r')
scatter(Temp_validation,TraPPE_validation,'c')
plot(Temp_exp,LDN_Funke)
plot(Temp_exp,LDN_DIPPR_exp,'m')
plot(Temp_exp,LDN_DIPPR_exp*1.01,'m--')
plot(Temp_exp,LDN_DIPPR_exp*0.99,'m--')
% % plot(Temp_max_min,rho_L_max,'r')
% % plot(Temp_max_min,rho_L_min,'r')
hold

% This is using the coarser grid
% epsilon = [95 95.5 96 96.5 97 98 99 99.5 100 102 102.5 103];
% sigma = [3.72 3.73 3.74 3.75 3.76 3.77 3.78];
% 
% All_data = dlmread('Coarse_All.txt');

% This is using the finer grid (missing 100.5 because of dump mistake)
epsilon = [96.5 97 97.5 98 98.5 99 99.5 100];
sigma = [3.72 3.725 3.73 3.735 3.74 3.745 3.75 3.755 3.76 3.765 3.77 3.775 3.78];

All_data = dlmread('Refined_All.txt');

% Rerefined grid
% epsilon = [97.5 97.75 98 98.25 98.5 98.75 99];
% sigma = [3.73 3.7325 3.735 3.7375 3.74 3.7425 3.745 3.7475 3.75 3.7525 3.755 3.7575 3.76 3.7625 3.765 3.7675 3.77];
% 
% All_data = dlmread('Rerefined_All.txt');

% Rererefined grid
% epsilon = [97.5 97.63 97.75 97.88 98 98.13 98.25 98.37 98.5 98.63 98.75 98.88]; % 99];
% sigma = [3.7375 3.7388 3.74 3.7413 3.7425 3.7438 3.745 3.7463 3.7475 3.7488 3.75 3.7513 3.7525 3.7538 3.755 3.7563 3.7575];

% 12 replicates region
% epsilon = [97.88 98 98.13 98.25 98.37 98.5 98.63 98.75 98.88];
% sigma = [3.745 3.7463 3.7475 3.7488 3.75 3.7513 3.7525 3.7538 3.755];

% 12 replicates region refined
% epsilon = [98.37 98.42 98.46 98.5 98.54 98.58 98.63 98.67 98.71 98.75];
% sigma = [3.7471 3.7475 3.7479 3.7483 3.7488 3.7492 3.7496 3.75 3.7504 3.7508 3.7513 3.7517 3.7521 3.7525];

liq_total = zeros(length(Temp),length(epsilon),length(sigma));
vap_total = liq_total;

liq_cat = zeros(1,length(epsilon),length(sigma));
vap_cat = liq_cat;
Temp_cat = 0;
TRC_cat = 0;
replicates = 0;

for run = 1:1 % This averages together multiple runs to get the true estimate from GEMC

%     All_data = dlmread(['Ethane_' int2str(run) '.txt']); % For the 12
    
%     All_data = dlmread(['Ethane_10' int2str(run) '.txt']); % For the 12 refined
    
    % I can think of some code to make these the same size grids or to
    % include the coarser grid in the analysis.
    if run < 4
        
    else
        
    end

for h = 1:length(epsilon) % Cycles through all the different epsilons

% eps_data = dlmread(sprintf('%s_All.txt',num2str(epsilon(h)))); % This is for when each different epsilon has a file

eps_data = All_data(:,(1+2*length(Temp)*(h-1)):(2*length(Temp)+2*length(Temp)*(h-1)));

[steps, temps] = size(eps_data);

temps = temps/2; % Since we only have half the columns (two per temp since one is liquid and one vapor)

steps = steps/length(sigma);

first_step = 1;
last_step = first_step + steps - 1;

for i = 1:length(sigma) % Cycles through all the different sigmas
    
    % Use this to investigate a single parameter set
     if epsilon(h) == 98.25 && sigma(i) == 3.7525
        
         pause(0)
         
     end
     
           
   for j = 1:temps
       
       swap = 0;
          
        for k = first_step:last_step
            
            dens_1 = eps_data(k,(2*j-1));
            dens_2 = eps_data(k,(2*j));
            
            if dens_1 >= dens_2  % This accounts for the possibility of identity swap between boxes
       
            liq_dens(k,j) = dens_1;
            vap_dens(k,j) = dens_2;
        
            else
                
            liq_dens(k,j) = dens_2;
            vap_dens(k,j) = dens_1;
            swap = k;
                
            end
            
        end
        
        equil_step = first_step + 2;
        
        repeat = true;
               
        while repeat == true
            
            if swap <= (equil_step)
        
                liq_ave(j,h,i) = mean(liq_dens(equil_step:last_step,j));
                vap_ave(j,h,i) = mean(vap_dens(equil_step:last_step,j));
        
                liq_std(j,h,i) = std(liq_dens(equil_step:last_step,j));
                vap_std(j,h,i) = std(vap_dens(equil_step:last_step,j));
        
        
                if liq_dens(equil_step,j) <= 0.9*liq_ave(j,h,i) ||  vap_dens(equil_step,j) >= 1.1*vap_ave(j,h,i)
        
            repeat = true;
            equil_step = equil_step + 1;
            
                else
            
            repeat = false;
        
                end
        
            else
                
                liq_swap = liq_dens(equil_step:swap-1,j);
                liq_swap = [liq_swap, liq_dens(swap+1:last_step,j)];
                vap_swap = vap_dens(equil_step:swap-1,j);
                vap_swap = [vap_swap, vap_dens(swap+1:last_step,j)];
                
                liq_ave(j,h,i) = mean(liq_swap);
                vap_ave(j,h,i) = mean(vap_swap);
        
                liq_std(j,h,i) = std(liq_swap);
                vap_std(j,h,i) = std(vap_swap);
        
        
                if liq_dens(equil_step,j) <= 0.9*liq_ave(j,h,i) ||  vap_dens(equil_step,j) >= 1.1*vap_ave(j,h,i)
        
            repeat = true;
            equil_step = equil_step + 1;
            
                else
            
            repeat = false;
        
                end
                
                
                
            end
        
        end
            
            
   end
    
   first_step = first_step + steps;
   last_step = first_step + steps - 1;
         
end     

end

% This is where I need to sum up the averages of replicates

liq_total = liq_total + liq_ave;
vap_total = vap_total + vap_ave;

liq_cat = [liq_cat; liq_ave];
vap_cat = [vap_cat; vap_ave];
Temp_cat = [Temp_cat, Temp];
TRC_cat = [TRC_cat, TRC];

replicates = replicates + 1;

end

% This just takes the average of the replicates
liq_ave_ave = liq_total / replicates;
vap_ave_ave = vap_total / replicates;

% This concatenates the replicates
liq_ave = liq_cat(2:end,:,:);
vap_ave = vap_cat(2:end,:,:);
Temp = Temp_cat(2:end);
TRC = TRC_cat(2:end);

% % % This shows the relationship between sigma, epsilon and rho_L
% 
% for j = 1:length(Temp_sim)
%     
%     X = ones(length(sigma),1);
%     X = [X, sigma'];
%     
% % figure
% 
% for h = 1:length(epsilon) % Span the different epsilons
%     
%     Y = liq_ave_ave(j,h,:);
%     Y = permute(Y,[3 2 1]);
%     
%     b_fit = (X'*X)\(X'*Y);
%     
%     Y_fit = X*b_fit;
%     
%     Y_sigma(j,h,:) = Y_fit;
%     
% %     hold
% %     scatter(sigma,Y)
% %     plot(sigma,Y_fit)
% %     title(Temp_sim(j))
% %     hold
%     
% end
% 
%     X = ones(length(epsilon),1);
%     X = [X, epsilon'];
%     
% % figure
% 
% for i = 1:length(sigma)
%     
%     Y = liq_ave_ave(j,:,i);
%     Y = permute(Y,[2 3 1]);
%     
%     Y_sig = Y_sigma(j,:,i);
%     Y_sig = permute(Y_sig,[2 3 1]);
%     
%     b_fit = (X'*X)\(X'*Y);
%     b_eps_sig = (X'*X)\(X'*Y_sig);
%     
%     Y_fit = X*b_fit;
%     Y_eps_sig = X*b_eps_sig;
%     
%     Y_epsilon(j,:,i) = Y_fit;
%     Y_epsilon_sigma(j,:,i) = Y_eps_sig;
%     
% %     hold
% %     scatter(epsilon,Y,'b')
% %     scatter(epsilon,Y_sig,'g')
% %     plot(epsilon,Y_fit,'b')
% %     plot(epsilon,Y_eps_sig,'g')
% %     title(Temp_sim(j))
% %     hold
% %     
% end
% 
% end

% global_minimum_ethane % If you want to do an advanced smoothing approach to find the true minimum and parameter confidence region, do this one

% reduced_ethane

valid = ones(length(epsilon),length(sigma),length(Temp));
valid95 = ones(length(epsilon),length(sigma));

   beta_range = linspace(0.326,0.326,1); % Seeing what the beta effect is
   
for h = 1:length(epsilon) % Cycles through all the different epsilons
       
   for i = 1:length(sigma) % Cycles through all the different sigmas
       
       % Use this to investigate a single parameter set
     if epsilon(h) == 98.46 && sigma(i) == 3.7483
        
         pause(0)
                 
     end
     
%      for j = 1:length(Temp_sim)
%    
%          for k = 1:replicates
%              
%             liq_sim(k) = liq_ave(j + (k-1)*n_sim,h,i);
%         
%          end
%          
%             std_sim(j,h,i) = std(liq_sim);
%             error_95_sim(j,h,i) =  2.22814 * std_sim(j,h,i) / sqrt(replicates);
%             per_error_95_sim(j,h,i) = error_95_sim(j,h,i) ./ liq_ave(j,h,i)' * 100;
% 
%      end
       
   for bb = 1:length(beta_range)
       
       beta = beta_range(bb);

   [TC(h,i), rhoc(h,i), A, b] = towhee_error_model(Temp,vap_ave(:,h,i)',liq_ave(:,h,i)',beta); % This uses my error model

%      [TC(h,i), rhoc(h,i), A,b] = towhee_rectilinear(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)',beta); % This is the traditional way
     
%     [rhoc_fit(h,i), TC_fit(h,i), A, b] = fit_rho_L(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta); % This just fits the rho_L data

%     [rhoc_fit(h,i), TC_fit(h,i), A, b, beta] = fit_rho_L_beta(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta); % This just fits the rho_L data

%     [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta); % This just fits the rho_L data but with the weighting model
    
      [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L(Temp_sim,liq_ave_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta); % This just fits the rho_L data but with the weighting model
      
%       [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L(Temp_sim,Y_epsilon_sigma(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta); % This fits the smoothed rho_L data with the weighting model

%     [b_Funke_fit] = weighted_fit_rho_L_Funke(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i));

%     TC_fit(h,i) = b_Funke_fit(2);
    
     liq_hat(:,h,i) = rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_exp) + b.*(TC_fit(h,i)-Temp_exp).^(beta);
     
%      liq_hat_Funke = @(T,b) b(1) .* exp(b(3) .* (1-T./b(2)).^(0.329) + b(4) .* (1-T./b(2)).^(4/6) + b(5) .* (1-T./b(2)).^(8/6) + b(6) .* (1-T./b(2)).^(19/6)); 

%      liq_hat(:,h,i) = liq_hat_Funke(Temp_exp,b_Funke_fit);
     
     liq_dev(:,h,i) = (rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_sim) + b.*(TC_fit(h,i)-Temp_sim).^(beta)) - liq_ave_ave(:,h,i)';
     
%      liq_dev(:,h,i) = liq_hat_Funke(Temp,b_Funke_fit) - liq_ave(:,h,i)';

     SSE_dev(h,i) = sum(liq_dev(:,h,i).^2);
     
     RMS_dev(h,i) = sqrt(SSE_dev(h,i)/n_sim);

     per_dev(:,h,i) = liq_dev(:,h,i) ./ liq_ave_ave(:,h,i) * 100;
     
     error_hat = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(Temp/DIPPR_TC));
     
     if TC_fit(h,i) > max(Temp_exp)
   
        SSE(h,i) = sum((DIPPR_exp'-liq_hat(:,h,i)).^2); % Non-weighted approach
%         SSE(h,i) = sum((DIPPR_exp'-liq_ave(:,h,i)).^2); % Non-weighted approach, using the average of the replicates (only valid when using TRC correlation)
%         SSE(h,i) = sum(((DIPPR_exp'-liq_hat(:,h,i))./error_hat).^2); % Weighted approach

%         SSE(h,i) = SSE(h,i) - SSE_dev(h,i); % Trying to account for the lack of fit to simulation
         
     else
         
        SSE(h,i) = 100000000;
         
     end
     
%         SSE_sim(h,i) = sum((TRC(1:6)'-liq_ave_ave(:,h,i)).^2); % Non-weighted approach, using average of replicate values
%         SSE(h,i) = sum(((TRC'-liq_ave(:,h,i))./error_hat').^2); % Weighted approach, using all the replicate values
%         SSE(h,i) = sum(((TRC(1:6)'-liq_ave_ave(:,h,i))./std_sim(:,h,i)).^2); % Weighted approach, using average of replicate values and the  simulation estimate of error (problematic because the estimates of error are very sporadic)
        SSE_sim(h,i) = sum((LDN_Funke_Sim'-liq_ave_ave(:,h,i)).^2); % Non-weighted approach, using all the replicate values

        
     if bb == 1
               
                SSE_opt(h,i) = SSE(h,i);
                beta_opt(h,i) = beta;
                TC_opt(h,i) = TC(h,i);
                liq_opt(:,h,i) = liq_hat(:,h,i);
                
     elseif SSE(h,i) < SSE_opt(h,i)
            
                SSE_opt(h,i) = SSE(h,i);
                beta_opt(h,i) = beta;
                TC_opt(h,i) = TC(h,i);
                liq_opt(:,h,i) = liq_hat(:,h,i);

     end
     
%      for ii = 1:length(Temp)
%          
%          % This does multiple z-statistic tests to determine how many are accepted
%          
%          z_stat = 1.115; % 0.991 (0.95^1/6)
%          t_stat = 1.175; % 3.196
%      
%         if sqrt((liq_ave(ii,h,i) - TRC(ii))^2) >= (t_stat / sqrt(replicates) * error_hat(ii))
%         
%             valid(h,i,ii) = 0;            
%             valid95(h,i) = 0; % This is 0 if any of the 6 temperatures is not accepted
%          
%         end
%         
%      end
     
   end
   
   end
   
end

% This tests each point to see if the regression lands within the 95%
% confidence region according to rho_L regression of experimental data, but
% remember this was using the entire density range

% z = 1;
% 
% for h = 1:length(epsilon)
%     
%     for i = 1:length(sigma)
%         
%         % This uses the average value of the replicates
%         upper_diff = max(rho_L_max-liq_ave(:,h,i)');
%         lower_diff = max(liq_ave(:,h,i)' - rho_L_min);
%         
%         % This uses the regression
% %         upper_diff = min(rho_L_max-liq_opt(:,h,i)');
% %         lower_diff = min(liq_opt(:,h,i)' - rho_L_min);
%         
%         if upper_diff > 0 && lower_diff > 0
%             
%             hold
%             plot(Temp_exp',liq_opt(:,h,i)','g')
%             scatter(Temp,liq_ave(:,h,i)','g')
%             hold
%             acceptable(:,z) = [epsilon(h) sigma(i)];
%             TC_acceptable(z) = TC_opt(h,i);
%             z=z+1;
%         
%         end
%         
%     end
%     
% end

% Trying to automate the contours
    p = 2; % Not sure what this value really should be, depends on the method
    n = length(Temp_exp);
%     n = length(Temp); % If the replicates are considered the data
    alpha = 0.95;

    [SSE_min, I] = min(SSE_opt);
    [SSE_min, J] = min(SSE_min);
    
    eps_opt = epsilon(I(J));
    sig_opt = sigma(J);
    opt = [eps_opt, sig_opt];
    
    RHS = SSE_min * (1 + (p/(n-p))*finv(alpha,p,n-p));
    
    sigma2 = SSE_min/(n-p);
    
    PDF = fcdf((SSE_opt - SSE_min)/(sigma2*p),p,n-p);
    
    RMS_opt = sqrt(SSE_opt/n);
    RMS_min = sqrt(SSE_min/n);
    RHS_RMS = sqrt(RHS/n);
    
    RMS_sim = sqrt(SSE_sim/n_sim);
    [RMS_sim_min, I] = min(RMS_sim);
    [RMS_sim_min, J] = min(RMS_sim_min);
    
    eps_opt_sim = epsilon(I(J));
    sig_opt_sim = sigma(J);
    opt_sim = [eps_opt_sim, sig_opt_sim];
    
    RHS_RMS_sim = RMS_sim_min * sqrt(1 + (p/(n-p))*finv(alpha,p,n-p));
    
    z = 1;
    zz = 1;
    
    RHS_alt = 2.09273363161464e-06;
       
for h = 1:length(epsilon)
    
    for i = 1:length(sigma)

%         if valid95(h,i) == 1 
        if SSE_opt(h,i) <= RHS_alt %SSE_opt(h,i) <= RHS %RMS_sim(h,i) <= RHS_RMS_sim
%             hold
%             plot(Temp_exp',liq_opt(:,h,i)','g')
%             scatter(Temp(1:6),liq_ave_ave(:,h,i)','g')
%             hold
            acceptable(:,z) = [epsilon(h) sigma(i)];
            TC_acceptable(z) = TC_opt(h,i);
            z=z+1;
%             hold
% %             scatter(Temp_sim,liq_dev(1:n_sim,h,i)','r')
%             scatter(Temp_sim,per_dev(1:n_sim,h,i)','r')
% %             plot(Temp_sim,std_sim(:,h,i))
% %             plot(Temp_sim,-std_sim(:,h,i))
%             plot(Temp_sim,per_error_95_sim(:,h,i),'g')
%             plot(Temp_sim,-per_error_95_sim(:,h,i),'g')
%             xlabel('Temperature (K)')
%             ylabel('Percent Deviation')
%             hold

        else
            
            not_acceptable(:,zz) = [epsilon(h) sigma(i)];
            zz = zz+1;
            
        end
        
%         scatter(Temp_sim,std_sim(:,h,i))
    
    end
    
end

error_hat_sim = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(Temp_sim/DIPPR_TC));

% hold
% plot(Temp_sim,error_hat_sim)
% hold

error_95 = 1.96 * error_hat_sim / sqrt(replicates);
per_error_95 = error_95 ./ liq_ave(1:n_sim,I(J),J)' * 100;

% hold
% plot(Temp_sim,per_error_95,'g')
% plot(Temp_sim,-per_error_95,'g')
% xlabel('Temperature (K)')
% ylabel('Percent Deviation')
% hold

%  sigma2 = 9.371*10^-10;
%  sigma2 = mean(error_hat)^2;
%  RHS_data = sigma2 * (n + p * (finv(0.95,p,n-p)-1)); % No 0.95^p because I want the true 95%

figure
contour(sigma,epsilon,RMS_opt)

% figure
% contour(sigma,epsilon,TC_opt)

% figure
% contour(sigma,epsilon,beta_opt)

% figure
% contour(sigma,epsilon,valid95)

% figure
% subplot(3,2,1)
% contour(sigma,epsilon,valid(:,:,1))
% subplot(3,2,2)
% contour(sigma,epsilon,valid(:,:,2))
% subplot(3,2,3)
% contour(sigma,epsilon,valid(:,:,3))
% subplot(3,2,4)
% contour(sigma,epsilon,valid(:,:,4))
% subplot(3,2,5)
% contour(sigma,epsilon,valid(:,:,5))
% subplot(3,2,6)
% contour(sigma,epsilon,valid(:,:,6))

figure
hold
scatter(acceptable(2,:),acceptable(1,:))
scatter(sig_opt,eps_opt,'r')
scatter(not_acceptable(2,:),not_acceptable(1,:),'g')
hold