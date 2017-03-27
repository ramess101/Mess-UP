clear

MW = 226.45;

Temp_sim = [510 625];
beta = 0.32;

n_sim = length(Temp_sim);

% From my error model, in case I need to estimate it
sig_L = [0.000836 0.000881 0.001146 0.001604 0.002719];
sig_v = [0.000356 0.000389 0.000661 0.001241 0.002906];

% Refined for optimal CH3
epsilon = [45.25 45.30 45.35 45.40 45.45 45.50];
sigma = [3.968 3.970 3.972 3.974 3.976 3.978 3.980 3.982];

liq_total = zeros(length(Temp_sim),length(epsilon),length(sigma));
vap_total = liq_total;

liq_cat = zeros(1,length(epsilon),length(sigma));
vap_cat = liq_cat;
Temp_cat = 0;
TRC_cat = 0;
replicates = 0;

for run = 1:5 % This averages together multiple runs to get the true estimate from GEMC
    
     All_data = dlmread(['C16_' int2str(run) '.txt']); % For the non-optimal runs

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

replicates = replicates + 1;

end

% This just takes the average of the replicates
liq_ave_ave = liq_total / replicates;
vap_ave_ave = vap_total / replicates;

% This concatenates the replicates
liq_ave = liq_cat(2:end,:,:);
vap_ave = vap_cat(2:end,:,:);
Temp_cat = Temp_cat(2:end);

% global_minimum_ethane % If you want to do an advanced smoothing approach to find the true minimum and parameter confidence region, do this one

for h = 1:length(epsilon) % Cycles through all the different epsilons
       
   for i = 1:length(sigma) % Cycles through all the different sigmas
      
   
   [TC(h,i), rhoc(h,i), A,b] = towhee_error_model(Temp_cat,vap_ave(:,h,i)',liq_ave(:,h,i)',beta); % This uses my error model
       
%    [rhoc_fit(h,i), TC_fit(h,i), A, b] = weighted_fit_rho_L_constrained(Temp_cat,liq_ave(:,h,i)',rhoc(h,i),TC(h,i),A,b,beta,max(Temp_sim)); % This just fits the rho_L data but keeps TC above the maximum experimental temperature
    
%    liq_hat(:,h,i) = rhoc_fit(h,i) + A.*(TC_fit(h,i)-Temp_exp) + b.*(TC_fit(h,i)-Temp_exp).^(beta); % This uses interpolation of just rhoL
   
   end
   
end

figure
contour(sigma,epsilon,TC)