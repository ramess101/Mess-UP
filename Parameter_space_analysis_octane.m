clear

% DIPPR = [0.624 0.574 0.505 0.473 0.425]; % DIPPR correlation
TRC = [0.620709305 0.569901486 0.507992855 0.469046647 0.418080225]; %TRC values
Temp = [390 440 490 515 540];
% DIPPR = [0.624 0.574 0.505 0.473];
% Temp = [390 440 490 515];
DIPPR_TC = 569;
beta = 0.32;

% hold
% plot(Temp,DIPPR)
% hold

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
% % Zoomed in grid
% % epsilon = [45 45.5 46 46.5 47];
% % sigma = [3.9 3.92 3.94 3.96 3.98 4.0];
% 
% All_data = dlmread('Octane_98.25_3.75_no_44.txt');

% Rezoomed in grid
% epsilon = [45.8 45.9 46 46.1 46.2];
% sigma = [3.945 3.948 3.95 3.952 3.955 3.957];

% Global search short runs
epsilon = [44.6 44.8 45 45.2 45.3 45.4 45.5 45.6 45.7 45.8 45.9 46 46.1 46.2 46.4 46.6 46.8 47];
sigma = [3.8 3.82 3.84 3.86 3.88 3.9 3.92 3.93 3.94 3.95 3.96 3.97 3.98 3.99 4.0 4.2 4.4];

All_data = dlmread('Octane_98_3.75_global_short.txt');


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
        
                liq_ave(j,i) = mean(liq_dens(equil_step:last_step,j));
                vap_ave(j,i) = mean(vap_dens(equil_step:last_step,j));
        
                liq_std(j,i) = std(liq_dens(equil_step:last_step,j));
                vap_std(j,i) = std(vap_dens(equil_step:last_step,j));
        

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
                liq_ave(j,i) = bin_liq(I_liq);
                liq_std(j,i) = sig_L(j); % Estimate the STD from the models I developed (at the exact TR)
                
                [n_vap, bin_vap] = hist(vap_dens(equil_step:last_step,j),20);
                [N_vap, I_vap] = max(n_vap);
                vap_ave(j,i) = bin_vap(I_vap);
                vap_std(j,i) = sig_v(j); % Estimate the STD from the models I developed (at the exact TR)
                
            end
            
                if liq_dens(equil_step,j) <= 0.9*liq_ave(j,i) ||  vap_dens(equil_step,j) >= 1.1*vap_ave(j,i)
        
            repeat = true;
            equil_step = equil_step + 1;
            
                else
            
            repeat = false;
        
                end
        
        end
            
            
   end
    
   first_step = first_step + steps;
   last_step = first_step + steps - 1;
   
   SSE(h,i) = sum((TRC'-liq_ave(:,i)).^2);
   
%    SSE(h,i) = sum(((DIPPR'-liq_ave(:,i))./(dLDN + liq_std(:,i))).^2); % This is a weighted SSE, however TraPPE used non-weighted
   
%    [TC(h,i), rhoc(h,i)] = towhee_error_model(Temp,vap_ave(:,i)',liq_ave(:,i)',beta); % This uses my error model

   [TC(h,i), rhoc(h,i)] = towhee_rectilinear(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)',beta); % This is the traditional way
   
   %    SETC(h,i) = ((DIPPR_TC-TC(h,i))/(dTC+TC_std(h,i))).^2; % This is a weigthed SSE
   
%    if TC(h,i) > 565 && TC(h,i) < 575
%    
%    [TC_low(h,i), TC_high(h,i), rhoc_low(h,i), rhoc_high(h,i)] = rigorous_statistics(Temp,vap_ave(:,i)',liq_ave(:,i)'); % This analyzes the uncertainty regions
%    
% %    [TC_low(h,i), TC_high(h,i), rhoc_low(h,i), rhoc_high(h,i)] = rigorous_statistics(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)'); % This analyzes the uncertainty regions
%       
%    end  
%     

%     if SSE(h,i) < 0.0003
%         hold
%         scatter(Temp,liq_ave(:,i)')
%         hold
%     end

end  
    

end

figure
contour(sigma,epsilon,SSE)
figure
contour(sigma,epsilon,TC)