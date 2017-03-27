clear

% DIPPR = [0.552729093 0.528132619 0.500108273 0.470700752 0.435279423 0.394205114]; % DIPPR correlation
TRC = [0.551945935 0.528130197 0.500847103 0.472019478 0.437027813 0.396053755]; % TRC correlation
Temp = [178 197 217 236 256 275];
DIPPR_TC = 305;
dLDN = 0.0012;
dTC = 0.04;

beta = 0.32;

hold
plot(Temp,TRC)
hold

% This is using the coarser grid
% epsilon = [95 95.5 96 96.5 97 98 99 99.5 100 102 102.5 103];
% sigma = [3.72 3.73 3.74 3.75 3.76 3.77 3.78];
% 
% All_data = dlmread('Coarse_All.txt');

% This is using the finer grid (missing 100.5 because of dump mistake)
% epsilon = [96.5 97 97.5 98 98.5 99 99.5 100];
% sigma = [3.72 3.725 3.73 3.735 3.74 3.745 3.75 3.755 3.76 3.765 3.77 3.775 3.78];

% Rerefined grid
epsilon = [97.5 97.75 98 98.25 98.5 98.75 99];
sigma = [3.73 3.7325 3.735 3.7375 3.74 3.7425 3.745 3.7475 3.75 3.7525 3.755 3.7575 3.76 3.7625 3.765 3.7675 3.77];

All_data = dlmread('Rerefined_All.txt');

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
        
        
                if liq_dens(equil_step,j) <= 0.9*liq_ave(j,i) ||  vap_dens(equil_step,j) >= 1.1*vap_ave(j,i)
        
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
                
                liq_ave(j,i) = mean(liq_swap);
                vap_ave(j,i) = mean(vap_swap);
        
                liq_std(j,i) = std(liq_swap);
                vap_std(j,i) = std(vap_swap);
        
        
                if liq_dens(equil_step,j) <= 0.9*liq_ave(j,i) ||  vap_dens(equil_step,j) >= 1.1*vap_ave(j,i)
        
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
   
   SSE(h,i) = sum((TRC'-liq_ave(:,i)).^2);
   
%    SSE(h,i) = sum(((DIPPR'-liq_ave(:,i))./(dLDN + liq_std(:,i))).^2); % This is a weighted SSE, however TraPPE used non-weighted
   
%    [TC(h,i), rhoc(h,i)] = towhee_error_model(Temp,vap_ave(:,i)',liq_ave(:,i)'); % This uses my error model

      [TC(h,i), rhoc(h,i)] = towhee_rectilinear(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)',beta); % This is the traditional way
   
   %    SETC(h,i) = ((DIPPR_TC-TC(h,i))/(dTC+TC_std(h,i))).^2; % This is a weigthed SSE
   
%    if TC(h,i) > 301.5 && TC(h,i) < 308
%    
% %    [TC_low(h,i), TC_high(h,i), rhoc_low(h,i), rhoc_high(h,i)] = rigorous_statistics(Temp,vap_ave(:,i)',liq_ave(:,i)'); % This analyzes the uncertainty regions
%    
%    [TC_low(h,i), TC_high(h,i), rhoc_low(h,i), rhoc_high(h,i)] = rigorous_statistics(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)'); % This analyzes the uncertainty regions
%       
%    end  
%     

    if SSE(h,i) < 0.00001 %&& h == 3 && i == 9
        hold
        scatter(Temp,liq_ave(:,i)')
        hold
%         epsilon(h)
%         sigma(i)
    end
        
end  
    

end

contour(sigma,epsilon,SSE)