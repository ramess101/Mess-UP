clear
% These are only the CH2 parameters
epsilon = [45.39479	45.37524	45.36170	45.40180	45.40782	45.33664	45.33965	45.45494	45.32160	45.39629	45.39880	45.40381	45.37674	45.41383	45.41033	45.37123	45.35419	45.34316	45.38426	45.38677	45.39128	45.30607	45.33263	45.37173	45.35619	45.38226	45.35368	45.39378	45.33815	45.44391	45.40080	45.36622	45.35268	45.37724	45.39078	45.37724	45.39028	45.35168	45.36070	45.39529	45.39278	45.39378	45.33514	45.37323	45.36271	45.38877	45.36221	45.37223	45.43338	45.36221	45.31208	45.40180	45.40882	45.37724	45.41333	45.32110	45.38977	45.37424	45.38276	45.40782	45.33764	45.36672	45.37073	45.40682	45.41283	45.37574	45.37925	45.37223	45.43840	45.40281	45.39529	45.38727	45.28150	45.33113	45.37975	45.40632	45.39278	45.35820	45.38526	45.38877	45.38576	45.36922	45.37223	45.38526	45.31609	45.35619	45.38376	45.37875	45.35970	45.35018	45.29704	45.38276	45.37424	45.39880	45.35820	45.40782	45.40030	45.35368	45.37023	45.34416];
sigma = [3.972699	3.976805	3.970353	3.973602	3.969947	3.972383	3.974549	3.974053	3.975947	3.972789	3.971256	3.969632	3.972068	3.970850	3.971256	3.973105	3.972564	3.974549	3.975677	3.970895	3.971301	3.979466	3.974729	3.972248	3.976083	3.973556	3.972338	3.973827	3.972880	3.967962	3.977752	3.976489	3.977887	3.970759	3.970534	3.970444	3.970805	3.975632	3.971211	3.971436	3.974414	3.970714	3.971165	3.972880	3.974729	3.969947	3.975090	3.972113	3.970850	3.975992	3.973872	3.971887	3.971752	3.972023	3.970714	3.976985	3.972338	3.973872	3.971346	3.970805	3.975677	3.975677	3.970805	3.972835	3.968143	3.969767	3.972699	3.974774	3.969722	3.974098	3.968323	3.971752	3.971301	3.975722	3.971842	3.971481	3.971165	3.972248	3.972880	3.973150	3.972203	3.974278	3.974549	3.970444	3.976038	3.973241	3.973150	3.976985	3.972519	3.977752	3.972970	3.972293	3.971075	3.969316	3.972383	3.969902	3.971571	3.970398	3.969767	3.973015];

MW = 506.97276;

Temp_exp = [473.15 523.15 573.15];
DIPPR_exp = [0.704793531 0.673006339 0.639901018];

% hold
% scatter(Temp_exp,DIPPR_exp)
% hold

% Temp_plot = linspace(min(Temp_exp),max(Temp_exp),1000);

alpha = 0.95;

% hold
% scatter(Temp_exp,DIPPR_exp,'k')
% hold

Temp = [650 775];
Temp_sim = Temp;
n_sim = length(Temp_sim);

liq_total = zeros(length(Temp),length(sigma));
vap_total = liq_total;
liq_std_total = liq_total;
vap_std_total = liq_total;

liq_cat = zeros(1,length(sigma));
vap_cat = liq_cat;
Temp_cat = 0;
replicates = 0;

for run = 1:5 % This averages together multiple runs to get the true estimate from GEMC
  
%     All_data = dlmread(['MC_C36_' int2str(run) '_eq.txt']); % For the Monte Carlo sampling sets
    
    All_data = dlmread(['MC_C36_Type_AB_' int2str(run) '.txt']); % For the Monte Carlo sampling sets
    
h = 1;

eps_data = All_data(:,(1+2*length(Temp)*(h-1)):(2*length(Temp)+2*length(Temp)*(h-1)));

[steps, temps] = size(eps_data);

temps = temps/2; % Since we only have half the columns (two per temp since one is liquid and one vapor)

steps = steps/length(sigma);

first_step = 1;
last_step = first_step + steps - 1;

for i = 1:length(sigma) % Cycles through all the different sigmas
    
    % Use this to investigate a single parameter set
     if epsilon(i) == 98.25 && sigma(i) == 3.7525
        
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
         
end     

% This is where I need to sum up the averages of replicates

liq_total = liq_total + liq_ave;
vap_total = vap_total + vap_ave;

liq_std_total = liq_std_total + liq_std;
vap_std_total = vap_std_total + vap_std;

liq_cat = [liq_cat; liq_ave];
vap_cat = [vap_cat; vap_ave];
Temp_cat = [Temp_cat, Temp];

replicates = replicates + 1;

end

% This just takes the average of the replicates
liq_ave_ave = liq_total / replicates;
vap_ave_ave = vap_total / replicates;

liq_std_ave = liq_std_total / replicates;
vap_std_ave = liq_std_total / replicates;

% This concatenates the replicates
liq_ave = liq_cat(2:end,:,:);
vap_ave = vap_cat(2:end,:,:);
Temp = Temp_cat(2:end);

% This shows the relationship between sigma, epsilon and rho_L
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
%     
% end
% 
% end

% global_minimum_ethane % If you want to do an advanced smoothing approach to find the true minimum and parameter confidence region, do this one

% reduced_ethane

valid = ones(length(sigma),length(Temp));
valid95 = ones(length(sigma));

beta = 0.32;
       
for i = 1:length(sigma) % Cycles through all the different sigmas
       
       % Use this to investigate a single parameter set
     if epsilon(i) == 98.54353 && sigma(i) == 3.749086
        
         pause(0)
                 
     end
%      
%      for j = 1:length(Temp_sim)
%    
%          for k = 1:replicates
%              
%             liq_sim(k) = liq_ave(j + (k-1)*n_sim,i);
%         
%          end
%          
%             std_sim(j,i) = std(liq_sim);
%             error_95_sim(j,i) =  2.22814 * std_sim(j,i) / sqrt(replicates);
%             per_error_95_sim(j,i) = error_95_sim(j,i) ./ liq_ave(j,i)' * 100;
% 
%      end
 
   [TC(i), rhoc(i), A, b] = towhee_error_model(Temp,vap_ave(:,i)',liq_ave(:,i)',beta); % This uses my error model

   [rhoc_fit, TC_fit, A, b] = weighted_fit_rho_L(Temp,liq_ave(:,i)',rhoc(i),TC(i),A,b,beta); % This fits the rho_L data with the weighting model
   
%    liq_hat(:,i) = rhoc_fit(i) + A.*(TC_fit(i)-Temp) + b.*(TC_fit(i)-Temp).^(beta);
%    
%    PC(i) = PC_from_Rackett_alt(TC(i),rhoc(i),liq_hat(1,i)',Temp(1)/TC(i),MW);

    liq_hat = rhoc_fit + A.*(TC_fit-Temp(1)) + b.*(TC_fit-Temp(1)).^(beta);
   
   PC(i) = PC_from_Rackett_alt(TC(i),rhoc(i),liq_hat,Temp(1)/TC(i),MW);
   
end

ZC = PC * MW ./ rhoc ./ TC ./ 8.314472;

for k = 1:length(Temp_sim)
    
    [rhoL_counts(k,:), rhoL_centers(k,:)] = hist(liq_ave_ave(k,:),20);

    [rhoL_lo_rig(k), rhoL_hi_rig(k)] = integrate_histogram(rhoL_counts(k,:),rhoL_centers(k,:),alpha);

    rhoL_uncertainty(k) = (rhoL_hi_rig(k) - rhoL_lo_rig(k)) / (rhoL_hi_rig(k) + rhoL_lo_rig(k)) * 100;

end

for k = 1:length(Temp_sim)
    
    [rhov_counts(k,:), rhov_centers(k,:)] = hist(vap_ave_ave(k,:),50);

    [rhov_lo_rig(k), rhov_hi_rig(k)] = integrate_histogram(rhov_counts(k,:),rhov_centers(k,:),alpha);

    rhov_uncertainty(k) = (rhov_hi_rig(k) - rhov_lo_rig(k)) / (rhov_hi_rig(k) + rhov_lo_rig(k)) * 100;

end

figure

for k =1:length(Temp_sim)
    
    subplot(2,2,k)
    hist(liq_ave_ave(k,:),20)

end


[TC_counts,TC_centers] = hist(TC,20);

[TC_lo_rig,TC_hi_rig] = integrate_histogram(TC_counts,TC_centers,alpha);

TC_uncertainty = (TC_hi_rig - TC_lo_rig) / (TC_hi_rig + TC_lo_rig) * 100;

figure
hist(TC,20)

[rhoc_counts,rhoc_centers] = hist(rhoc,20);

[rhoc_lo_rig,rhoc_hi_rig] = integrate_histogram(rhoc_counts,rhoc_centers,alpha);

rhoc_uncertainty = (rhoc_hi_rig - rhoc_lo_rig) / (rhoc_hi_rig + rhoc_lo_rig) * 100;

figure
hist(rhoc,20)

[PC_counts,PC_centers] = hist(PC,20);

[PC_lo_rig,PC_hi_rig] = integrate_histogram(PC_counts,PC_centers,alpha);

PC_uncertainty = (PC_hi_rig - PC_lo_rig) / (PC_hi_rig + PC_lo_rig) * 100;

figure
hist(PC,20)

[ZC_counts,ZC_centers] = hist(ZC,20);

[ZC_lo_rig,ZC_hi_rig] = integrate_histogram(ZC_counts,ZC_centers,alpha);

ZC_uncertainty = (ZC_hi_rig - ZC_lo_rig) / (ZC_hi_rig + ZC_lo_rig) * 100;

figure
hist(ZC,20)

rhoL_centers(1,1)
rhoL_centers(1,2) - rhoL_centers(1,1)
rhoL_centers(2,1)
rhoL_centers(2,2) - rhoL_centers(2,1)
TC_centers(1)
(TC_centers(20)-TC_centers(1))/20
rhoc_centers(1)
(rhoc_centers(20)-rhoc_centers(1))/20
PC_centers(1)
(PC_centers(20)-PC_centers(1))/20
ZC_centers(1)
(ZC_centers(20)-ZC_centers(1))/20

criticals=[TC',rhoc',PC',ZC'];

1.96*std(liq_ave_ave(1,:))
1.96*std(liq_ave_ave(2,:))