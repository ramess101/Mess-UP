clear

MW = 114.23;

DIPPR_exp = [0.616993928 0.611990719 0.607992721 0.598991513 0.587488701 0.589990306 0.57999531 0.5769911 0.568995104 0.560896302 0.557994897 0.54598948 0.531288269 0.533995485 0.530991275 0.511995072 0.508990862 0.497796467 0.494997869 0.478994453 0.472997456 0.455394841 0.45999825 0.437997837 0.410994215]; 
Temp_exp = [393.15 398.15 403.15 413.15 423.15 423.15 433.15 433.15 443.15 448.15 453.15 463.15 473.15 473.15 473.15 483.15 493.15 498.15 503.15 513.15 513.15 523.15 523.15 533.15 543.15]; 
Temp_TRC = [393.15 403.15 413.15 423.15 433.15 443.15 453.15 463.15 473.15 483.15 493.15 503.15 513.15 523.15 533.15 543.15];
TRC = [0.616993928	0.607615766	0.597997725	0.588094112	0.5777907	0.567190293	0.555995898	0.544287475	0.531893681	0.518688864	0.504490259	0.489297865	0.474493849	0.455988829	0.434993627	0.407990005];
Temp_acc = [393.15 398.15 403.15 413.15 423.15 433.15 448.15 473.15 498.15 513.15 523.15 543.15];
LDN_acc = [0.616993928 0.611990719 0.607992721 0.598991513 0.587488701 0.5769911 0.560896302 0.530991275 0.497796467 0.472997456 0.455394841 0.410994215];

Temp_sim = [390 440 490 515 540];
TraPPE = [0.624 0.574 0.505 0.473 0.425];

DIPPR_TC = 569;
beta = 0.32;

n_sim = length(Temp_sim);

% DIPPR_exp = DIPPR_exp(Temp_exp<535);
% Temp_exp = Temp_exp(Temp_exp<535);
% DIPPR_exp = DIPPR_exp(Temp_exp~=483.15);
% Temp_exp = Temp_exp(Temp_exp~=483.15);

% DIPPR_exp = TRC;
% Temp_exp = Temp_TRC;

% LDN_acc = [LDN_acc DIPPR_exp(length(DIPPR_exp)-1)];
% Temp_acc = [Temp_acc Temp_exp(length(Temp_exp)-1)];

DIPPR_exp = LDN_acc;
Temp_exp = Temp_acc;

DIPPR_exp = DIPPR_exp(Temp_exp<540);
Temp_exp = Temp_exp(Temp_exp<540);

hold
scatter(Temp_exp,DIPPR_exp)
hold

Parameters_C = [0.213240223 568.6487323 2.254331272 -1.285838279 0.584843942 -0.127688134]; % Fit to all the rho_L data

Funke_fit = Parameters_C;

Funke_hat = @(b,T) b(1) * exp(b(3)*(1-T/b(2)).^0.329 + b(4)*(1-T/b(2)).^(4/6)+b(5)*(1-T/b(2)).^(8/6) + b(6)*(1-T/b(2)).^(19/6));

LDN_Funke = @(T) Funke_hat(Funke_fit,T);

Temp_plot = linspace(min(Temp_exp),max(Temp_exp),1000);

LDN_Funke_exp = LDN_Funke(Temp_exp);
LDN_Funke_Sim = LDN_Funke(Temp_sim);

LDN_Funke_plot = LDN_Funke(Temp_plot);

SSE_Funke = sum((DIPPR_exp'-LDN_Funke_exp').^2); % Non-weighted approach

 p = 6;
 n = length(Temp_exp);
 
 sigma2_exp = SSE_Funke/(n-p);
 
 alpha = 0.95;

 RHS_Funke = SSE_Funke * (1 + (p/(n-p))*finv(alpha,p,n-p));
 
hold
scatter(Temp_exp,DIPPR_exp,'k')
% scatter(Temp_TRC,TRC,'r')
scatter(Temp_sim,TraPPE,'c')
plot(Temp_plot,LDN_Funke_plot,'b')
hold

Temp = Temp_sim;

% 100 randomly chosen locations with 8 replicates
epsilon_100 = [45.40884	45.37957	45.38616	45.37711	45.38838	45.35758	45.35996	45.34644	45.37518	45.37309	45.40273	45.35886	45.34655	45.40764	45.41076	45.38745	45.36350	45.38083	45.41063	45.40379	45.41809	45.33494	45.40000	45.41825	45.33462	45.41090	45.33233	45.34106	45.36191	45.43153	45.40245	45.43357	45.44398	45.33853	45.44771	45.40194	45.41017	45.41000	45.39398	45.34691	45.31803	45.34011	45.30529	45.32000	45.42571	45.37976	45.35507	45.44470	45.42523	45.34525	45.30643	45.44310	45.32341	45.42733	45.31954	45.41879	45.41837	45.43607	45.43248	45.41599	45.32500	45.43849	45.33600	45.34410	45.46526	45.46165	45.34519	45.29049	45.29342	45.36015	45.39076	45.33500	45.30967	45.34635	45.44215	45.48498	45.45997	45.37082	45.45517	45.46296	45.28905	45.36026	45.45241	45.41415	45.38994	45.41983	45.28029	45.41215	45.31000	45.29992	45.39155	45.28848	45.45845	45.47209	45.46198	45.38482	45.31029	45.28229	45.36676	45.44062];
sigma_100 = [3.972541	3.974855	3.973226	3.973240	3.974336	3.975679	3.977363	3.977071	3.976392	3.977034	3.973629	3.976570	3.975697	3.974492	3.973249	3.974983	3.975280	3.975974	3.973754	3.972930	3.974085	3.975507	3.977000	3.971767	3.977946	3.972021	3.977604	3.978225	3.974127	3.973269	3.971336	3.970519	3.970995	3.975328	3.971354	3.970861	3.971395	3.976000	3.971507	3.973799	3.978380	3.979455	3.978238	3.979500	3.971188	3.972613	3.978801	3.972032	3.975293	3.979073	3.977904	3.973738	3.977019	3.974975	3.976685	3.975156	3.976013	3.973636	3.970850	3.971042	3.979000	3.975800	3.973999	3.974121	3.969485	3.972277	3.980495	3.977055	3.980106	3.979788	3.977849	3.980000	3.975864	3.979963	3.974118	3.970191	3.969132	3.978782	3.972732	3.971327	3.980479	3.972262	3.968980	3.969873	3.970506	3.970308	3.980156	3.969299	3.979000	3.981346	3.970187	3.979473	3.973584	3.969331	3.970886	3.971184	3.980181	3.980611	3.971189	3.969111];

% The next 53 samples
epsilon_53 = [45.29073	45.30091	45.30504	45.31000	45.31024	45.31689	45.32000	45.32056	45.32194	45.32225	45.33022	45.33335	45.33860	45.34000	45.34324	45.34430	45.34499	45.35252	45.35267	45.35274	45.35347	45.35359	45.35524	45.35633	45.35991	45.36498	45.36735	45.36958	45.36970	45.37108	45.37109	45.37500	45.37543	45.37691	45.38000	45.39041	45.39079	45.39243	45.39262	45.40044	45.40525	45.41000	45.41567	45.41640	45.42500	45.42501	45.42646	45.42920	45.42953	45.44000	45.44027	45.44554	45.44609];
sigma_53 =[3.977785	3.975229	3.973450	3.974500	3.976516	3.971212	3.975000	3.973391	3.970486	3.971587	3.973818	3.971073	3.971643	3.972500	3.969242	3.976143	3.968128	3.971619	3.974684	3.973447	3.972909	3.970130	3.967745	3.968289	3.966879	3.969984	3.972463	3.969054	3.968175	3.973979	3.969812	3.966000	3.974642	3.972142	3.967000	3.973991	3.968126	3.969228	3.965046	3.967212	3.966481	3.968000	3.967047	3.966306	3.969000	3.965926	3.965445	3.967455	3.967980	3.967000	3.965003	3.968872	3.964125];

eps_ext = [45.37957	45.37711	45.35758	45.35996	45.37309	45.35886	45.34655	45.40764	45.41076	45.38745	45.36350	45.38083	45.41809	45.33494	45.33462	45.41090	45.34106	45.36191	45.40245	45.43357	45.33853	45.40194	45.41017	45.42571	45.32341	45.31954	45.43248	45.41599	45.30967	45.36026	45.41983	45.39155	45.38482	45.31000	45.31024	45.32000	45.32056	45.33022	45.33335	45.35274	45.35347	45.35359	45.36498	45.36958	45.39079	45.39243	45.40525	45.41640	45.42500	45.42920	45.42953	45.44554];
sig_ext = [3.974855	3.973240	3.975679	3.977363	3.977034	3.976570	3.975697	3.974492	3.973249	3.974983	3.975280	3.975974	3.974085	3.975507	3.977946	3.972021	3.978225	3.974127	3.971336	3.970519	3.975328	3.970861	3.971395	3.971188	3.977019	3.976685	3.970850	3.971042	3.975864	3.972262	3.970308	3.970187	3.971184	3.974500	3.976516	3.975000	3.973391	3.973818	3.971073	3.973447	3.972909	3.970130	3.969984	3.969054	3.968126	3.969228	3.966481	3.966306	3.969000	3.967455	3.967980	3.968872];

epsilon = [epsilon_100 epsilon_53];
sigma = [sigma_100 sigma_53];

epsilon_surrogate = 45.37697;
sigma_surrogate = 3.972273;

liq_total = zeros(length(Temp),length(sigma));
vap_total = liq_total;

liq_cat = zeros(1,length(sigma));
vap_cat = liq_cat;
Temp_cat = 0;
TRC_cat = 0;
replicates = 0;

for run = 1:9 % This averages together multiple runs to get the true estimate from GEMC
    
%     if run > 3 % Run 4 has some weird issue with the data analysis
%         
%         run = run + 1;
%                                
%     end
      
    All_data = dlmread(['Octane_ellipse_' int2str(run) '.txt']); % For the 9 ellipse
    
h = 1;

eps_data = All_data(:,(1+2*length(Temp)*(h-1)):(2*length(Temp)+2*length(Temp)*(h-1)));

octane_extended_data_analysis

[steps, temps] = size(eps_data);

temps = temps/2; % Since we only have half the columns (two per temp since one is liquid and one vapor)

steps = steps/length(sigma);

first_step = 1;
last_step = first_step + steps - 1;

for i = 1:length(sigma) % Cycles through all the different sigmas
    
    % Use this to investigate a single parameter set
     if i == 28
        
         pause(0)
         
     end
     
           
   for j = 1:temps
       
       swap = 0;
          
        for k = first_step:last_step
            
            dens_1 = eps_data(k,(2*j-1));
            dens_2 = eps_data(k,(2*j));
            
            if dens_1 >= dens_2  % This accounts for the possibility of identity swap between boxes
       
                if swap == 0
                
                    liq_dens(k,j) = dens_1;
                    vap_dens(k,j) = dens_2;
                                        
                end
                
                swap = 0;
        
            else
                
                if swap == 1
                
                    liq_dens(k,j) = dens_2;
                    vap_dens(k,j) = dens_1;
                                
                end
                
                swap = 1;
                
            end
            
        end
        
        equil_step = first_step + 2;
        
        repeat = true;
               
        while repeat == true
            
            if swap <= (equil_step)
                
               liq_sample = liq_dens(equil_step:last_step,j);
               liq_sample = liq_sample(liq_sample~=0);
                
               vap_sample = vap_dens(equil_step:last_step,j);
               vap_sample = vap_sample(vap_sample~=0);
                
               liq_ave(j,i) = mean(liq_sample);
               vap_ave(j,i) = mean(vap_sample);
                
               liq_std(j,i) = std(liq_sample);
               vap_std(j,i) = std(vap_sample);
%         
%                  liq_ave(j,i) = mean(liq_dens(equil_step:last_step,j));
%                  vap_ave(j,i) = mean(vap_dens(equil_step:last_step,j));
%          
%                  liq_std(j,i) = std(liq_dens(equil_step:last_step,j));
%                  vap_std(j,i) = std(vap_dens(equil_step:last_step,j));
                
                outlier = true;
                
                while outlier == true % Make sure to eliminate outliers
                    
                    outlier = false;
                    
                    liq_old = length(liq_sample);
                    vap_old = length(vap_sample);
                    
                    liq_sample = liq_sample(liq_sample <= liq_ave(j,i) + 3*liq_std(j,i));
                    liq_sample = liq_sample(liq_sample >= liq_ave(j,i) - 3*liq_std(j,i));
                    
                    vap_sample = vap_sample(vap_sample <= vap_ave(j,i) + 3*vap_std(j,i));
                    vap_sample = vap_sample(vap_sample >= vap_ave(j,i) - 3*vap_std(j,i));
                    
                    liq_new = length(liq_sample);
                    vap_new = length(vap_sample);
                    
                    if abs(liq_old - liq_new) > 0 || abs(vap_old - vap_new) > 0
                        
                        outlier = true;
                        
                    end
                    
                end
                    
%                     t = 1;
%                     tt = 1;
%                    
%                     for kk = equil_step:last_step
%                        
%                         if liq_dens(kk,j) <= liq_ave(j,i) + 3*liq_std(j,i) && liq_dens(kk,j) >= liq_ave(j,i) - 3*liq_std(j,i)
%                             
%                             liq_sample(t) = liq_dens(kk,j);
%                             
%                             t = t + 1;
%                             
%                         else
%                             
%                             outlier = true;
%                                                      
%                         end
%                         
%                         if vap_dens(kk,j) <= vap_ave(j,i) + 3*vap_std(j,i) && vap_dens(kk,j) >= vap_ave(j,i) - 3*vap_std(j,i)
%                             
%                             vap_sample(tt) = vap_dens(kk,j);
%                             
%                             tt = tt + 1;
%                             
%                         else
%                             
%                             outlier = true;
%                                                      
%                         end
%                         
%                     end
                    
                    liq_ave(j,i) = mean(liq_sample);
                    vap_ave(j,i) = mean(vap_sample);
                    
                    liq_std(j,i) = std(liq_sample);
                    vap_std(j,i) = std(vap_sample);

%                 clear liq_sample vap_sample
%                     
%                 end
        
        
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

   beta_range = linspace(0.326,0.326,1); % Seeing what the beta effect is
       
for i = 1:length(sigma) % Cycles through all the different sigmas
       
       % Use this to investigate a single parameter set
     if epsilon(i) == 98.46 && sigma(i) == 3.7483
        
         pause(0)
                 
     end
     
     for j = 1:length(Temp_sim)
   
         for k = 1:replicates
             
            liq_sim(k) = liq_ave(j + (k-1)*n_sim,i);
        
         end
         
            std_sim(j,i) = std(liq_sim);
            error_95_sim(j,i) =  2.22814 * std_sim(j,i) / sqrt(replicates);
            per_error_95_sim(j,i) = error_95_sim(j,i) ./ liq_ave(j,i)' * 100;

     end
       
   for bb = 1:length(beta_range)
       
       beta = beta_range(bb);

   [TC(i), rhoc(i), A, b] = towhee_error_model(Temp,vap_ave(:,i)',liq_ave(:,i)',beta); % This uses my error model

%      [TC(i), rhoc(i), A,b] = towhee_rectilinear(Temp,vap_ave(:,i)',liq_ave(:,i)',vap_std(:,i)',liq_std(:,i)',beta); % This is the traditional way
     
%     [rhoc_fit(i), TC_fit(i), A, b] = fit_rho_L(Temp,liq_ave(:,i)',rhoc(i),TC(i),A,b,beta); % This just fits the rho_L data

%     [rhoc_fit(i), TC_fit(i), A, b, beta] = fit_rho_L_beta(Temp,liq_ave(:,i)',rhoc(i),TC(i),A,b,beta); % This just fits the rho_L data

%     [rhoc_fit(i), TC_fit(i), A, b] = weighted_fit_rho_L(Temp,liq_ave(:,i)',rhoc(i),TC(i),A,b,beta); % This just fits the rho_L data but with the weighting model
    
      [rhoc_fit(i), TC_fit(i), A, b] = weighted_fit_rho_L(Temp_sim,liq_ave_ave(:,i)',rhoc(i),TC(i),A,b,beta); % This just fits the rho_L data but with the weighting model
      
%       [rhoc_fit(i), TC_fit(i), A, b] = weighted_fit_rho_L(Temp_sim,Y_epsilon_sigma(:,i)',rhoc(i),TC(i),A,b,beta); % This fits the smoothed rho_L data with the weighting model

%     [b_Funke_fit] = weighted_fit_rho_L_Funke(Temp,liq_ave(:,i)',rhoc(i),TC(i));

%     TC_fit(i) = b_Funke_fit(2);
    
     liq_hat_exp(:,i) = rhoc_fit(i) + A.*(TC_fit(i)-Temp_exp) + b.*(TC_fit(i)-Temp_exp).^(beta);
     
     liq_hat_sim(:,i) = rhoc_fit(i) + A.*(TC_fit(i)-Temp_sim) + b.*(TC_fit(i)-Temp_sim).^(beta);
     
     liq_hat_plot(:,i) = rhoc_fit(i) + A.*(TC_fit(i)-Temp_plot) + b.*(TC_fit(i)-Temp_plot).^(beta);
     
%      liq_hat_Funke = @(T,b) b(1) .* exp(b(3) .* (1-T./b(2)).^(0.329) + b(4) .* (1-T./b(2)).^(4/6) + b(5) .* (1-T./b(2)).^(8/6) + b(6) .* (1-T./b(2)).^(19/6)); 

%      liq_hat(:,i) = liq_hat_Funke(Temp_exp,b_Funke_fit);
     
     liq_dev(:,i) = (rhoc_fit(i) + A.*(TC_fit(i)-Temp_sim) + b.*(TC_fit(i)-Temp_sim).^(beta)) - liq_ave_ave(:,i)';
     
%      liq_dev(:,i) = liq_hat_Funke(Temp,b_Funke_fit) - liq_ave(:,i)';

     SSE_dev(i) = sum(liq_dev(:,i).^2);
          
     RMS_dev(i) = sqrt(SSE_dev(i)/n_sim);

     per_dev(:,i) = liq_dev(:,i) ./ liq_ave_ave(:,i) * 100;
     
     error_hat = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(Temp/DIPPR_TC));
     
     if TC_fit(i) > max(Temp_exp)
   
        SSE(i) = sum((DIPPR_exp'-liq_hat_exp(:,i)).^2); % Non-weighted approach
%         SSE(i) = sum((DIPPR_exp'-liq_ave(:,i)).^2); % Non-weighted approach, using the average of the replicates (only valid when using TRC correlation)
%         SSE(i) = sum(((DIPPR_exp'-liq_hat(:,i))./error_hat).^2); % Weighted approach

%         SSE(i) = SSE(i) - SSE_dev(i); % Trying to account for the lack of fit to simulation
         
     else
         
        SSE(i) = 100000000;
         
     end
     
%         SSE_sim(i) = sum((TRC(1:6)'-liq_ave_ave(:,i)).^2); % Non-weighted approach, using average of replicate values
%         SSE(i) = sum(((TRC'-liq_ave(:,i))./error_hat').^2); % Weighted approach, using all the replicate values
%         SSE(i) = sum(((TRC(1:6)'-liq_ave_ave(:,i))./std_sim(:,i)).^2); % Weighted approach, using average of replicate values and the  simulation estimate of error (problematic because the estimates of error are very sporadic)
        SSE_sim(i) = sum((LDN_Funke_Sim'-liq_ave_ave(:,i)).^2); % Non-weighted approach, using all the replicate values

        
     if bb == 1
               
                SSE_opt(i) = SSE(i);
                beta_opt(i) = beta;
                TC_opt(i) = TC(i);
                liq_opt(:,i) = liq_hat_exp(:,i);
                liq_opt_plot(:,i) = liq_hat_plot(:,i);
                
     elseif SSE(i) < SSE_opt(i)
            
                SSE_opt(i) = SSE(i);
                beta_opt(i) = beta;
                TC_opt(i) = TC(i);
                liq_opt(:,i) = liq_hat_exp(:,i);
                liq_opt_plot(:,i) = liq_hat_plot(:,i);

     end
     
%      for ii = 1:length(Temp)
%          
%          % This does multiple z-statistic tests to determine how many are accepted
%          
%          z_stat = 1.115; % 0.991 (0.95^1/6)
%          t_stat = 1.175; % 3.196
%      
%         if sqrt((liq_ave(ii,i) - TRC(ii))^2) >= (t_stat / sqrt(replicates) * error_hat(ii))
%         
%             valid(i,ii) = 0;            
%             valid95(i) = 0; % This is 0 if any of the 6 temperatures is not accepted
%          
%         end
%         
%      end
     
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
%         upper_diff = max(rho_L_max-liq_ave(:,i)');
%         lower_diff = max(liq_ave(:,i)' - rho_L_min);
%         
%         % This uses the regression
% %         upper_diff = min(rho_L_max-liq_opt(:,i)');
% %         lower_diff = min(liq_opt(:,i)' - rho_L_min);
%         
%         if upper_diff > 0 && lower_diff > 0
%             
%             hold
%             plot(Temp_exp',liq_opt(:,i)','g')
%             scatter(Temp,liq_ave(:,i)','g')
%             hold
%             acceptable(:,z) = [epsilon(i) sigma(i)];
%             TC_acceptable(z) = TC_opt(i);
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

    sigma2_sim = SSE_dev/(n_sim-p);

    alpha = 0.95;
    
    SSE_opt_alt = SSE_opt - SSE_dev * n / (n_sim - p); % Trying to see if the model fit is part of the problem
    
    [SSE_sort, J] = sort(SSE_opt);
    
    eps_sort = epsilon(J);
    sig_sort = sigma(J);
    
    SSE_opt_alt_sort = SSE_opt_alt(J);

    [SSE_min, I] = min(SSE_opt);
    
    eps_opt = epsilon(I);
    sig_opt = sigma(I);
    opt = [eps_opt, sig_opt];
    
%     SSE_min = sort(SSE_sim); % This uses SSE_sim, not SSE_opt
%     SSE_min = SSE_min(1); % Just in case there is an anomaly
    
    % Sort(*) is the value that corresponds best to the surrogate model
    SSE_min = SSE_sort(7); % Just in case you don't want the best fit
    eps_opt = eps_sort(7);
    sig_opt = sig_sort(7);

%     SSE_min = SSE_opt_alt_sort(1); % Alternative attempt
    
    RHS = SSE_min * (1 + (p/(n-p))*finv(alpha,p,n-p));
    
    RHS_sort = SSE_sort * (1 + (p/(n-p))*finv(alpha,p,n-p));
    
%     RHS = 1.59*10^-6;
    
    RMS_opt = sqrt(SSE_opt/n);
    RMS_min = sqrt(SSE_min/n);
    RHS_RMS = sqrt(RHS/n);
    
    RMS_sim = sqrt(SSE_sim/n_sim);
    [RMS_sim_min, I] = min(RMS_sim);
        
    eps_opt_sim = epsilon(I);
    sig_opt_sim = sigma(I);
    opt_sim = [eps_opt_sim, sig_opt_sim];
    
    RMS_sim_min = sort(RMS_sim);
    RMS_sim_min = RMS_sim_min(2);
    
    RHS_RMS_sim = RMS_sim_min * sqrt(1 + (p/(n-p))*finv(alpha,p,n-p));
    
%     sigma2 = sigma2_exp + max(sigma2_sim);
%     sigma2 = sigma2_exp + (mean(mean(std_sim)))^2;
%     sigma2 = sigma2_exp + 0.0014^2; % This value is the average of the error model for rho_L
    sigma2 = sigma2_exp + SSE_min/(n-p);
    
    RHS_alt = sigma2 * (n + p * (finv(alpha,p,n-p)-1));
    
    RHS_alt_2 = SSE_min + sigma2*p*finv(alpha,p,n-p);
    
    RHS_alt_3 = SSE_Funke + sigma2*p*finv(alpha,p,n-p);
    
    z = 1;
    zz = 1;
           
for i = 1:length(sigma)

%         if valid95(i) == 1 
        if  SSE_opt(h,i) <= RHS % SSE_opt(h,i) <= SSE_sort(30) %SSE_opt(h,i) <= RHS %RMS_sim(h,i) <= RHS_RMS_sim % SSE_opt(h,i) <= RHS_alt
            hold
            plot(Temp_plot',liq_opt_plot(:,i)','g')
            scatter(Temp(1:5),liq_ave_ave(:,i)','g')
            hold
            acceptable(:,z) = [epsilon(i) sigma(i)];
            TC_acceptable(z) = TC_opt(i);
            z=z+1;
%             hold
% %             scatter(Temp_sim,liq_dev(1:n_sim,i)','r')
%             scatter(Temp_sim,per_dev(1:n_sim,i)','r')
% %             plot(Temp_sim,std_sim(:,i))
% %             plot(Temp_sim,-std_sim(:,i))
%             plot(Temp_sim,per_error_95_sim(:,i),'g')
%             plot(Temp_sim,-per_error_95_sim(:,i),'g')
%             xlabel('Temperature (K)')
%             ylabel('Percent Deviation')
%             hold

        else
            
            not_acceptable(:,zz) = [epsilon(i) sigma(i)];
            TC_not_acceptable(zz) = TC_opt(i);
            hold
            plot(Temp_plot',liq_opt_plot(:,i)','r')
            scatter(Temp(1:5),liq_ave_ave(:,i)','r')
            hold
            
            zz = zz+1;
            
        end
        
%         scatter(Temp_sim,std_sim(:,i))
    
end

error_hat_sim = (8.27*10^(-4))+ (837*10^(-14))*exp(20.98*(Temp_sim/DIPPR_TC));

% hold
% plot(Temp_sim,error_hat_sim)
% hold

error_95 = 1.96 * error_hat_sim / sqrt(replicates);
per_error_95 = error_95 ./ liq_ave(1:n_sim,I)' * 100;

% hold
% plot(Temp_sim,per_error_95,'g')
% plot(Temp_sim,-per_error_95,'g')
% xlabel('Temperature (K)')
% ylabel('Percent Deviation')
% hold

%  sigma2 = 9.371*10^-10;
%  sigma2 = mean(error_hat)^2;
%  RHS_data = sigma2 * (n + p * (finv(0.95,p,n-p)-1)); % No 0.95^p because I want the true 95%

% figure
% contour(sigma,epsilon,RMS_opt)

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
scatter(acceptable(2,:),acceptable(1,:),'b')
scatter(not_acceptable(2,:),not_acceptable(1,:),'r')
scatter(sig_opt,eps_opt,'k')
scatter(sigma_surrogate,epsilon_surrogate,'p')
hold

figure
hold
plot(sort(SSE_opt),'b')
plot(RHS_alt*SSE_opt./SSE_opt,'r')
plot(RHS_alt_2*SSE_opt./SSE_opt,'g')
plot(RHS_alt_3*SSE_opt./SSE_opt,'c')
plot(RHS*SSE_opt./SSE_opt,'k')
plot(RHS_sort,'m')
hold

% figure
% hold
% scatter(sigma,epsilon,'k')
% scatter(sigma_surrogate,epsilon_surrogate,'p')
% for i = 1:length(sig_sort)
%     scatter(sig_sort(1:i),eps_sort(1:i),'r')
%     pause(1)
%     i
% end
