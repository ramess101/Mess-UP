
Mw_ethane = 30.07; 

k = 1;

for j = 1:length(Temp_sim)
    
    for h = 1:length(epsilon)
        
        for i = 1:length(sigma)
            
            rhoL_reduced(j,h,i) = liq_ave_ave(j,h,i) * (sigma(i)^3) / Mw_ethane * 0.6022; % This should be in reduced units with sigma in Angstroms and liq_ave_ave in gm/mL
            
            rhov_reduced(j,h,i) = vap_ave_ave(j,h,i) * (sigma(i)^3) / Mw_ethane * 0.6022;
            
            T_reduced(j,h,i) = Temp_sim(j) / epsilon(h);
            
                if epsilon(h) >= 98.42 && epsilon(h) <= 98.63 && sigma(i) >= 3.7479 && sigma(i) <= 3.7508
            
                    rhoL_reduced_all(k) = rhoL_reduced(j,h,i);
            
                    rhov_reduced_all(k) = rhov_reduced(j,h,i);
            
                    T_reduced_all(k) = T_reduced(j,h,i);
                        
                    k = k+1;
            
                end
            
        end
        
        
    end
    
end

hold
scatter(rhoL_reduced_all,T_reduced_all)
scatter(rhov_reduced_all,T_reduced_all)
hold