clear

MW = 114.23;

DIPPR_exp = [0.616993928 0.611990719 0.607992721 0.598991513 0.587488701 0.589990306 0.57999531 0.5769911 0.568995104 0.560896302 0.557994897 0.54598948 0.531288269 0.533995485 0.530991275 0.511995072 0.508990862 0.497796467 0.494997869 0.478994453 0.472997456 0.455394841 0.45999825 0.437997837 0.410994215]; 
Temp_exp = [393.15 398.15 403.15 413.15 423.15 423.15 433.15 433.15 443.15 448.15 453.15 463.15 473.15 473.15 473.15 483.15 493.15 498.15 503.15 513.15 513.15 523.15 523.15 533.15 543.15]; 
Temp_TRC = [393.15 403.15 413.15 423.15 433.15 443.15 453.15 463.15 473.15 483.15 493.15 503.15 513.15 523.15 533.15 543.15];
TRC = [0.616993928	0.607615766	0.597997725	0.588094112	0.5777907	0.567190293	0.555995898	0.544287475	0.531893681	0.518688864	0.504490259	0.489297865	0.474493849	0.455988829	0.434993627	0.407990005];
Temp_acc = [393.15 398.15 403.15 413.15 423.15 433.15 448.15 473.15 498.15 513.15 523.15 543.15];
LDN_acc = [0.616993928 0.611990719 0.607992721 0.598991513 0.587488701 0.5769911 0.560896302 0.530991275 0.497796467 0.472997456 0.455394841 0.410994215];

Temp = [390 440 490 515 540];
DIPPR_TC = 569;
beta = 0.32;

% DIPPR_exp = DIPPR_exp(Temp_exp<535);
% Temp_exp = Temp_exp(Temp_exp<535);
DIPPR_exp = DIPPR_exp(Temp_exp~=483.15);
Temp_exp = Temp_exp(Temp_exp~=483.15);

% DIPPR_exp = TRC;
% Temp_exp = Temp_TRC;

DIPPR_exp = LDN_acc;
Temp_exp = Temp_acc;

hold
scatter(Temp_exp,DIPPR_exp)
hold

% From my error model, in case I need to estimate it
sig_L = [0.000836 0.000881 0.001146 0.001604 0.002719];
sig_v = [0.000356 0.000389 0.000661 0.001241 0.002906];

% This is using the coarser grid
% epsilon = [95 95.5 96 96.5 97 98 99 99.5 100 102 102.5 103];
% sigma = [3.72 3.73 3.74 3.75 3.76 3.77 3.78];
% 
% All_data = dlmread('Coarse_All.txt');

% This is using the finer grid (missing 100.5 because of dump mistake)
% epsilon = [96.5 97 97.5 98 98.5 99 99.5 100];
% sigma = [3.72 3.725 3.73 3.735 3.74 3.745 3.75 3.755 3.76 3.765 3.77 3.775 3.78];

% Rerefined grid
% epsilon = [45 46 47 48 49];
% sigma = [3.94 3.945 3.95 3.955 3.96];
% 
% All_data = dlmread('Octane_98.25_3.75_no_44.txt');

% % Zoomed in grid
% epsilon = [45 45.5 46 46.5 47];
% sigma = [3.9 3.92 3.94 3.96 3.98 4.0];
% 
% All_data = dlmread('Octane_98.25_3.75_zoomed.txt');

% Rezoomed in grid
% epsilon = [45.8 45.9 46 46.1 46.2];
% sigma = [3.945 3.948 3.95 3.952 3.955 3.957];
% 
% All_data = dlmread('Octane_98.25_3.75_rezoomed.txt');

% Global search short runs
% epsilon = [44.6 44.8 45 45.2 45.3 45.4 45.5 45.6 45.7 45.8 45.9 46 46.1 46.2 46.4 46.6 46.8 47];
% sigma = [3.8 3.82 3.84 3.86 3.88 3.9 3.92 3.93 3.94 3.95 3.96 3.97 3.98 3.99 4.0 4.2 4.4];
% 
% All_data = dlmread('Octane_98_3.75_global_short.txt');

% Refined short runs
epsilon = [45.35 45.45 45.55 45.65 45.75 45.85 45.95 46.05];
sigma = [3.93 3.935 3.94 3.945 3.95 3.955 3.96 3.965 3.97 3.975 3.98 3.985 3.99 3.995 4.0];

All_data = dlmread('Octane_98_3.75_1.txt');



for h = 1:length(epsilon) % Cycles through all the different epsilons

% eps_data = dlmread(sprintf('%s_All.txt',num2str(epsilon(h)))); % This is for when each different epsilon has a file

eps_data = All_data(:,(1+2*length(Temp)*(h-1)):(2*length(Temp)+2*length(Temp)*(h-1)));

[steps, temps] = size(eps_data);

temps = temps/2; % Since we only have half the columns (two per temp since one is liquid and one vapor)

steps = steps/length(sigma);

first_step = 1;
last_step = first_step + steps - 1;

for i = 1:length(sigma) % Cycles through all the different sigmas
           
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
        

            else
                
                % My original attempt to deal with swapping, however this
                % performs poorly if a few outliers persist
                
%                 liq_swap = liq_dens(equil_step:swap-1,j);
%                 liq_swap = [liq_swap; liq_dens(swap+1:last_step,j)];
%                 vap_swap = vap_dens(equil_step:swap-1,j);
%                 vap_swap = [vap_swap; vap_dens(swap+1:last_step,j)];
%                 
%                 liq_ave(j,i) = mean(liq_swap);
%                 vap_ave(j,i) = mean(vap_swap);
%         
%                 liq_std(j,i) = std(liq_swap);
%                 vap_std(j,i) = std(vap_swap);     

                % I decided to adopt the histogram binning method
                
                [n_liq,bin_liq] = hist(liq_dens(equil_step:last_step,j),20);
                [N_liq,I_liq] = max(n_liq);
                liq_ave(j,h,i) = bin_liq(I_liq);
                liq_std(j,h,i) = sig_L(j); % Estimate the STD from the models I developed (at the exact TR)
                
                [n_vap, bin_vap] = hist(vap_dens(equil_step:last_step,j),20);
                [N_vap, I_vap] = max(n_vap);
                vap_ave(j,h,i) = bin_vap(I_vap);
                vap_std(j,h,i) = sig_v(j); % Estimate the STD from the models I developed (at the exact TR)
                
            end
            
                if liq_dens(equil_step,j) <= 0.9*liq_ave(j,h,i) ||  vap_dens(equil_step,j) >= 1.1*vap_ave(j,h,i)
        
            repeat = true;
            equil_step = equil_step + 1;
            
                else
            
            repeat = false;
        
                end
        
        end
            
            
   end
    
   first_step = first_step + steps;
   last_step = first_step + steps - 1;
   
   [TC(h,i), rhoc(h,i), A,b] = towhee_error_model(Temp,vap_ave(:,h,i)',liq_ave(:,h,i)',beta); % This uses my error model

%      [TC(h,i), rhoc(h,i), A,b] = towhee_rectilinear(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)',beta); % This is the traditional way

%      [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta); % This just fits the rho_L data
        
%      [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L_constrained(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta,max(Temp_exp)); % This just fits the rho_L data but keeps TC above the maximum experimental temperature

     [rhoc_fit(h,i), TC_fit(h,i), A, b] = fit_rho_L_constrained(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta,max(Temp_exp)); % This just fits the rho_L data but keeps TC above the maximum experimental temperature without weighting the data
     
%      liq_hat(:,h,i) = rhoc(h,i) + A.*(TC(h,i)-Temp_exp) + b.*(TC(h,i)-Temp_exp).^(beta); % This uses the rhor rhos
     
     liq_hat(:,h,i) = rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_exp) + b.*(TC_fit(h,i)-Temp_exp).^(beta); % This uses interpolation of just rhoL
   
%      if h >= 5 && h <= 13 && i >= 8 && i <=14
%      
%      figure
%      hold
%      plot(Temp_exp',liq_hat(:,h,i)')
%      scatter(Temp,liq_ave(:,h,i)')
%      hold
%      
%      end
     
     SSE(h,i) = sum((DIPPR_exp'-liq_hat(:,h,i)).^2);
   
%    if TC(h,i) > 565 && TC(h,i) < 575
%    
%    [TC_low(h,i), TC_high(h,i), rhoc_low(h,i), rhoc_high(h,i)] = rigorous_statistics(Temp,vap_ave(:,i)',liq_ave(:,i)'); % This analyzes the uncertainty regions
%    
% %    [TC_low(h,i), TC_high(h,i), rhoc_low(h,i), rhoc_high(h,i)] = rigorous_statistics(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)'); % This analyzes the uncertainty regions
%       
%    end  
%     

%     if SSE(h,i) < 0.0001976
%         hold
%         plot(Temp_exp',liq_hat(:,h,i)')
%         scatter(Temp,liq_ave(:,h,i)')
%         hold
%         acceptable(:,z) = [epsilon(h) sigma(i)];
%         z=z+1;
%     end

end  
    

end

% Trying to automate the contours
    p = 2;
    n = length(Temp_exp);
    alpha = 0.95;

    [SSE_min, I] = min(SSE);
    [SSE_min, J] = min(SSE_min);
    
    eps_opt = epsilon(I(J));
    sig_opt = sigma(J);
    opt = [eps_opt, sig_opt];
    
    RHS = SSE_min * (1 + (p/(n-p))*finv(alpha,p,n-p));
    
    RMS = sqrt(SSE/n);
    RMS_min = sqrt(SSE_min/n);
    RMS_RHS = sqrt(RHS/n);
    
    z = 1;
    
for h = 1:length(epsilon)
    
    for i = 1:length(sigma)

        if SSE(h,i) <= RHS
            hold
            plot(Temp_exp',liq_hat(:,h,i)')
            scatter(Temp,liq_ave(:,h,i)')
            hold
            acceptable(:,z) = [epsilon(h) sigma(i)];
            z=z+1;
        end
    
    end
    
end

figure
contour(sigma,epsilon,RMS)
figure
contour(sigma,epsilon,TC)
figure
hold
scatter(acceptable(2,:),acceptable(1,:))
scatter(sig_opt,eps_opt,'r')

Temp_sim = Temp;
n_sim = length(Temp_sim);
liq_ave_ave=liq_ave;
vap_ave_ave=vap_ave;