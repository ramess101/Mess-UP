clear

Temp_exp = [88.71 90.37 91 94.26 95 99.82 100 105.37 110 110.93 116.48 120 122.04 127.59 130 133.15 140 150 160 170 175 185 195 200 210 220 230 240 250 260 265 270 275 280 283.2 283.2 285 288.19 288.19 290 293.18 293.18 295 295.67 295.67 298 298.17 298.17 300 300.66 300.66 301 302 302.16 302.16 303 303.16 304.15 304.65 305.15 303.5 304 304.4 304.7 304.9 305 305.1 305.2 305.25 305.275]; 
DIPPR_exp = [0.6519176 0.65152669 0.650847108 0.64611409 0.646468916 0.64040079 0.640987155 0.63450707 0.630014612 0.62882384 0.62290005 0.618981929 0.61700633 0.6110224 0.60786505 0.60509861 0.59661887 0.585207305 0.573594271 0.561722635 0.555675558 0.543319795 0.530558087 0.52399982 0.51046832 0.496284301 0.48130042 0.465318215 0.448052021 0.429086872 0.418766848 0.407731158 0.39581141 0.382745995 0.37392045 0.37355961 0.368150017 0.35810363 0.35768265 0.351340887 0.33952037 0.33891897 0.33098049 0.32851475 0.32773293 0.31582521 0.31579514 0.31489304 0.30352658 0.300104614 0.29880559 0.296322409 0.28803632 0.288004446 0.286684373 0.278054884 0.276346307 0.262565226 0.252503804 0.235375932 0.272034 0.264902 0.257857 0.251177 0.245477 0.241946 0.237573 0.231703 0.227394 0.224471];    
Temp = [178 197 217 236 256 275];
DIPPR_TC = 305;
dLDN = 0.0012;
dTC = 0.04;
beta = 0.326;

% DIPPR_exp = DIPPR_exp(Temp_exp<=275);
% Temp_exp = Temp_exp(Temp_exp<=275);
% DIPPR_exp = DIPPR_exp(Temp_exp>=175);
% Temp_exp = Temp_exp(Temp_exp>=175);

hold
scatter(Temp_exp,DIPPR_exp)
hold

z = 1;

% This is using the coarser grid
% epsilon = [95 95.5 96 96.5 97 98 99 99.5 100 102 102.5 103];
% sigma = [3.72 3.73 3.74 3.75 3.76 3.77 3.78];
% 
% All_data = dlmread('Coarse_All.txt');

% This is using the finer grid (missing 100.5 because of dump mistake)
% epsilon = [96.5 97 97.5 98 98.5 99 99.5 100];
% sigma = [3.72 3.725 3.73 3.735 3.74 3.745 3.75 3.755 3.76 3.765 3.77 3.775 3.78];

% Rerefined grid
% epsilon = [97.5 97.75 98 98.25 98.5 98.75 99];
% sigma = [3.73 3.7325 3.735 3.7375 3.74 3.7425 3.745 3.7475 3.75 3.7525 3.755 3.7575 3.76 3.7625 3.765 3.7675 3.77];
% 
% All_data = dlmread('Rerefined_All.txt');

% Rererefined grid
epsilon = [97.5 97.63 97.75 97.88 98 98.13 98.25 98.37 98.5 98.63 98.75 98.88 99];
sigma = [3.7375 3.7388 3.74 3.7413 3.7425 3.7438 3.745 3.7463 3.7475 3.7488 3.75 3.7513 3.7525 3.7538 3.755 3.7563 3.7575];

All_data = dlmread('Rererefined_all_1.txt');

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
   
   beta_range = linspace(0.32,0.32,1); % Seeing what the beta effect is
   
   for bb = 1:length(beta_range)
       
       beta = beta_range(bb);
     
   [TC(h,i), rhoc(h,i), A, b] = towhee_error_model(Temp,vap_ave(:,h,i)',liq_ave(:,h,i)',beta); % This uses my error model

%      [TC(h,i), rhoc(h,i), A,b] = towhee_rectilinear(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)',beta); % This is the traditional way
     
     liq_hat(:,h,i) = rhoc(h,i) + A.*(TC(h,i)-Temp_exp) + b.*(TC(h,i)-Temp_exp).^(beta);
     
     error_hat(:,i) = (8.27*10^(-4))+ (8.37*10^(-14))*exp(20.98*(Temp_exp/DIPPR_TC));
     
     if TC(h,i) > max(Temp_exp)
   
        SSE(h,i) = sum((DIPPR_exp'-liq_hat(:,h,i)).^2); % Non-weighted approach
%         SSE(h,i) = sum(((DIPPR_exp'-liq_hat(:,i))./error_hat(:,i)).^2); % Weighted approach
        
     else
         
        SSE(h,i) = 100000000;
         
     end
     
     if bb == 1
               
                SSE_opt(h,i) = SSE(h,i);
                beta_opt(h,i) = beta;
                TC_opt(h,i) = TC(h,i);
                liq_opt(:,i) = liq_hat(:,i);
                
     elseif SSE(h,i) < SSE_opt(h,i)
            
                SSE_opt(h,i) = SSE(h,i);
                beta_opt(h,i) = beta;
                TC_opt(h,i) = TC(h,i);
                liq_opt(:,i) = liq_hat(:,i);

     end
     
   end
   
%    if TC(h,i) > 301.5 && TC(h,i) < 308
%    
% %    [TC_low(h,i), TC_high(h,i), rhoc_low(h,i), rhoc_high(h,i)] = rigorous_statistics(Temp,vap_ave(:,i)',liq_ave(:,i)'); % This analyzes the uncertainty regions
%    
%    [TC_low(h,i), TC_high(h,i), rhoc_low(h,i), rhoc_high(h,i)] = rigorous_statistics(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)'); % This analyzes the uncertainty regions
%       
%    end  
%    

%         if SSE_opt(h,i) <= (1.5653*10^-4)
%             hold
%             plot(Temp_exp',liq_opt(:,i)')
%             scatter(Temp,liq_ave(:,i)')
%             hold
%             acceptable(:,z) = [epsilon(h) sigma(i)];
%             z=z+1;
%         end
      
end     

end

% Trying to automate the contours
    p = 2;
    n = length(Temp_exp);
    alpha = 0.95;

    SSE_min = min(min(SSE_opt));
    RHS = SSE_min * (1 + (p/(n-p))*finv(alpha,p,n-p));
    
%     z = 1;
%     
% for h = 1:length(epsilon)
%     
%     for i = 1:length(sigma)
% 
%         if SSE(h,i) <= RHS
%             hold
%             plot(Temp_exp',liq_hat(:,i)')
%             hold
%             acceptable(:,z) = [epsilon(h) sigma(i)];
%             z=z+1;
%         end
%     
%     end
%     
% end

% figure
% contour(sigma,epsilon,SSE)