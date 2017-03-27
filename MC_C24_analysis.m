clear
% These are only the CH2 parameters
epsilon = [45.37724	45.35168	45.37223	45.39679	45.35268	45.39679	45.34717	45.37875	45.37825	45.32561	45.43338	45.36421	45.42486	45.38977	45.38927	45.41584	45.34967	45.35970	45.38376	45.39579	45.38175	45.38727	45.40381	45.41734	45.41083	45.40130	45.37624	45.35218	45.34416	45.37424	45.32912	45.32060	45.34266	45.36471	45.42185	45.37474	45.47148	45.29454	45.35368	45.30456	45.36471	45.40932	45.33764	45.41383	45.39479	45.37624	45.34617	45.40231	45.34216	45.36271	45.32712	45.36120	45.37474	45.37875	45.37273	45.37925	45.36070	45.37474	45.41434	45.39679	45.37774	45.42135	45.37323	45.39429	45.32160	45.37774	45.37624	45.40732	45.34216	45.37774	45.35318	45.36421	45.31459	45.38276	45.35018	45.36622	45.37875	45.32712	45.36070	45.35068	45.36922	45.32010	45.37825	45.33063	45.39980	45.34216	45.34967	45.38376	45.35769	45.33414	45.40782	45.39028	45.40982	45.34416	45.39830	45.39880	45.39880	45.37875	45.30256	45.44692];
sigma = [3.969992	3.973466	3.974684	3.971211	3.975451	3.973195	3.974684	3.973556	3.970669	3.974910	3.968684	3.970534	3.973150	3.970759	3.971617	3.974684	3.975947	3.975316	3.970985	3.975857	3.974729	3.974098	3.968459	3.970308	3.972158	3.970895	3.973647	3.972429	3.971481	3.972203	3.973376	3.975406	3.970398	3.974368	3.974188	3.971120	3.964398	3.979602	3.972564	3.977842	3.972474	3.968504	3.973602	3.971436	3.972383	3.972158	3.976038	3.971707	3.972113	3.973962	3.972519	3.971436	3.972744	3.973015	3.968774	3.973376	3.972699	3.970714	3.968008	3.974549	3.972338	3.969451	3.972654	3.972744	3.974820	3.974504	3.973195	3.972023	3.974053	3.972203	3.974729	3.972744	3.975722	3.972654	3.972023	3.972023	3.977120	3.974188	3.974278	3.975541	3.973060	3.975677	3.975496	3.976534	3.966925	3.973376	3.971842	3.971707	3.973331	3.975947	3.970579	3.973692	3.970850	3.976985	3.970308	3.973511	3.971436	3.972293	3.978970	3.972248];

MW = 338.65;

% hold
% scatter(Temp_exp,DIPPR_exp)
% hold

% Temp_plot = linspace(min(Temp_exp),max(Temp_exp),1000);

alpha = 0.95;

% hold
% scatter(Temp_exp,DIPPR_exp,'k')
% hold

Temp = [600 700];
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
  
%     All_data = dlmread(['MC_C24_' int2str(run) '.txt']); % For the Monte Carlo sampling sets
    
    All_data = dlmread(['MC_C24_Type_AB_' int2str(run) '.txt']); % For the Monte Carlo sampling sets Type AB
    
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

ZC = PC * MW / 8.314472 ./ TC ./ rhoc;

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