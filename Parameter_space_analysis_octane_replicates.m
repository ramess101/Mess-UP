clear

MW = 114.23;

DIPPR_exp = [0.616993928 0.611990719 0.607992721 0.598991513 0.587488701 0.589990306 0.57999531 0.5769911 0.568995104 0.560896302 0.557994897 0.54598948 0.531288269 0.533995485 0.530991275 0.511995072 0.508990862 0.497796467 0.494997869 0.478994453 0.472997456 0.455394841 0.45999825 0.437997837 0.410994215]; 
Temp_exp = [393.15 398.15 403.15 413.15 423.15 423.15 433.15 433.15 443.15 448.15 453.15 463.15 473.15 473.15 473.15 483.15 493.15 498.15 503.15 513.15 513.15 523.15 523.15 533.15 543.15]; 
Temp_TRC = [393.15 403.15 413.15 423.15 433.15 443.15 453.15 463.15 473.15 483.15 493.15 503.15 513.15 523.15 533.15 543.15];
TRC = [0.616993928	0.607615766	0.597997725	0.588094112	0.5777907	0.567190293	0.555995898	0.544287475	0.531893681	0.518688864	0.504490259	0.489297865	0.474493849	0.455988829	0.434993627	0.407990005];
Temp_acc = [393.15 398.15 403.15 413.15 423.15 433.15 448.15 473.15 498.15 513.15 523.15 543.15];
LDN_acc = [0.616993928 0.611990719 0.607992721 0.598991513 0.587488701 0.5769911 0.560896302 0.530991275 0.497796467 0.472997456 0.455394841 0.410994215];

Temp_sim = [390 440 490 515 540];
DIPPR_TC = 569;
beta = 0.32;

n_sim = length(Temp_sim);

% DIPPR_exp = DIPPR_exp(Temp_exp<535);
% Temp_exp = Temp_exp(Temp_exp<535);
% DIPPR_exp = DIPPR_exp(Temp_exp~=483.15);
% Temp_exp = Temp_exp(Temp_exp~=483.15);

% DIPPR_exp = TRC;
% Temp_exp = Temp_TRC;

DIPPR_exp = LDN_acc;
Temp_exp = Temp_acc;

DIPPR_exp = DIPPR_exp(Temp_exp<540);
Temp_exp = Temp_exp(Temp_exp<540);

% Parameters_DIPPR = [0.5266 0.25693 568.7 0.28571]; % This is the actual DIPPR correlation
Parameters_DIPPR = [0.53278965 0.2599501 568.7 0.2795282]; % This is my fit of Rackett to data keeping DIPPR rhoc and TC
% Parameters_DIPPR = [0.53658038 0.2607781 568.7 0.2809242]; % This is my fit of Rackett to data keeping DIPPR TC but actual VC (slightly different)
 
LDN_DIPPR_exp = MW / 1000 * Parameters_DIPPR(1) ./ (Parameters_DIPPR(2) .^ (1+(1-Temp_exp/Parameters_DIPPR(3)).^Parameters_DIPPR(4)));
LDN_DIPPR_sim = MW / 1000 * Parameters_DIPPR(1) ./ (Parameters_DIPPR(2) .^ (1+(1-Temp_sim/Parameters_DIPPR(3)).^Parameters_DIPPR(4)));

hold
scatter(Temp_exp,DIPPR_exp)
plot(Temp_exp,LDN_DIPPR_exp)
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
epsilon = [45 45.5 46 46.5 47];
sigma = [3.9 3.92 3.94 3.96 3.98 4.0];

All_data = dlmread('Octane_98.25_3.75_zoomed.txt');

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

% % Refined short runs
% epsilon = [45.35 45.45 45.55 45.65 45.75 45.85 45.95 46.05];
% sigma = [3.93 3.935 3.94 3.945 3.95 3.955 3.96 3.965 3.97 3.975 3.98 3.985 3.99 3.995 4.0];

% Refined for non-optimal CH3
% epsilon = [45.20 45.30 45.40 45.50 45.60];
% sigma = [3.965 3.970 3.975 3.980 3.985];

liq_total = zeros(length(Temp_sim),length(epsilon),length(sigma));
vap_total = liq_total;

liq_cat = zeros(1,length(epsilon),length(sigma));
vap_cat = liq_cat;
Temp_cat = 0;
TRC_cat = 0;
replicates = 0;

for run = 1:1 % This averages together multiple runs to get the true estimate from GEMC
    
%     All_data = dlmread(['Octane_98_3.75_' int2str(run) '.txt']); % For the 8 zoomed in runs at TraPPE
%     All_data = dlmread(['Octane_98.4966_3.7491_' int2str(run) '.txt']); % For the 9 zoomed in runs
%         All_data = dlmread(['Octane_98.395_3.748_' int2str(run) '.txt']); % For the non-optimal runs

for h = 1:length(epsilon) % Cycles through all the different epsilons

% eps_data = dlmread(sprintf('%s_All.txt',num2str(epsilon(h)))); % This is for when each different epsilon has a file

eps_data = All_data(:,(1+2*length(Temp_sim)*(h-1)):(2*length(Temp_sim)+2*length(Temp_sim)*(h-1)));

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
   
end

end

% This is where I need to sum up the averages of replicates

liq_total = liq_total + liq_ave;
vap_total = vap_total + vap_ave;

liq_cat = [liq_cat; liq_ave];
vap_cat = [vap_cat; vap_ave];
Temp_cat = [Temp_cat, Temp_sim];
TRC_cat = [TRC_cat, TRC];

replicates = replicates + 1;

end

% This just takes the average of the replicates
liq_ave_ave = liq_total / replicates;
vap_ave_ave = vap_total / replicates;

% This concatenates the replicates
liq_ave = liq_cat(2:end,:,:);
vap_ave = vap_cat(2:end,:,:);
Temp_cat = Temp_cat(2:end);
TRC = TRC_cat(2:end);

% global_minimum_ethane % If you want to do an advanced smoothing approach to find the true minimum and parameter confidence region, do this one

for h = 1:length(epsilon) % Cycles through all the different epsilons
       
   for i = 1:length(sigma) % Cycles through all the different sigmas

             % Use this to investigate a single parameter set
             
             if epsilon(h) == 98.46 && sigma(i) == 3.7483
                 
                 pause(0)
                 
             end
     
             % Problem child from TraPPE
%      if h == 8 && i == 2 
%         
%          pause(0)
%          liq_ave_ave(5,h,i) = (liq_ave(40,h,i) + liq_ave(35,h,i) + liq_ave(30,h,i) + liq_ave(20,h,i)+liq_ave(15,h,i)+liq_ave(10,h,i)+liq_ave(5,h,i))/7;
%          liq_ave(25,h,i) = liq_ave_ave(5,h,i);
%          vap_ave_ave(5,h,i) = (vap_ave(40,h,i) + vap_ave(35,h,i) + vap_ave(30,h,i) + vap_ave(20,h,i)+vap_ave(15,h,i)+vap_ave(10,h,i)+vap_ave(5,h,i))/7;
%          vap_ave(25,h,i) = vap_ave_ave(5,h,i);
%                  
%      end
     
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
       
   
   [TC(h,i), rhoc(h,i), A,b] = towhee_error_model(Temp_cat,vap_ave(:,h,i)',liq_ave(:,h,i)',beta); % This uses my error model

%      [TC(h,i), rhoc(h,i), A,b] = towhee_rectilinear(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)',beta); % This is the traditional way

%      [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta); % This just fits the rho_L data
        
     [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L_constrained(Temp_cat,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta,max(Temp_exp)); % This just fits the rho_L data but keeps TC above the maximum experimental temperature

%      [rhoc_fit(h,i), TC_fit(h,i), A, b] = fit_rho_L_constrained(Temp,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta,max(Temp_exp)); % This just fits the rho_L data but keeps TC above the maximum experimental temperature without weighting the data
     
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

     liq_dev(:,h,i) = (rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_sim) + b.*(TC_fit(h,i)-Temp_sim).^(beta)) - liq_ave_ave(:,h,i)';

     SSE_dev(h,i) = sum(liq_dev(:,h,i).^2);
     
     RMS_dev(h,i) = sqrt(SSE_dev(h,i)/n_sim);

     per_dev(:,h,i) = liq_dev(:,h,i) ./ liq_ave_ave(:,h,i) * 100;
     
     error_hat = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(Temp_sim/DIPPR_TC));
     
     SSE(h,i) = sum((DIPPR_exp'-liq_hat(:,h,i)).^2);
     
     SSE_sim(h,i) = SSE_dev(h,i);   

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

SSE_opt = SSE;
TC_opt = TC;
liq_opt = liq_hat;


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
        if SSE_opt(h,i) <= RHS %SSE_opt(h,i) <= RHS %RMS_sim(h,i) <= RHS_RMS_sim SSE_opt(h,i) <= RHS_alt
            hold
            plot(Temp_exp',liq_opt(:,h,i)','g')
            scatter(Temp_sim(1:5),liq_ave_ave(:,h,i)','g')
            hold
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