
Mw_octane = 114.23; 

k = 1;

for j = 1:length(Temp)
    
    for h = 1:length(epsilon)
        
        for i = 1:length(sigma)
            
            rhoL_reduced(j,h,i) = liq_ave(j,h,i) * (sigma(i)^3) / Mw_octane * 0.6022; % This should be in reduced units with sigma in Angstroms and liq_ave_ave in gm/mL
            
            rhov_reduced(j,h,i) = vap_ave(j,h,i) * (sigma(i)^3) / Mw_octane * 0.6022;
            
            T_reduced(j,h,i) = Temp(j) / epsilon(h);
            
            rhoL_reduced_all(k) = rhoL_reduced(j,h,i);
            
            rhov_reduced_all(k) = rhov_reduced(j,h,i);
            
            T_reduced_all(k) = T_reduced(j,h,i);
                        
            k = k+1;
            
        end
        
        
    end
    
end

hold
scatter(rhoL_reduced_all,T_reduced_all)
scatter(rhov_reduced_all,T_reduced_all)
hold